/* cmio.c
 * SRE, Thu Aug  3 11:53:34 2000 [St. Louis]
 * SVN $Id$
 * 
 * Input/output of covariance models.
 * 
 *****************************************************************
 * @LICENSE@
 ***************************************************************** 
 */ 

#include "config.h"

#include <stdio.h>
#include <string.h>
#include "squid.h"
#include "ssi.h"

#include "structs.h"
#include "funcs.h"

/* Magic numbers identifying binary formats.
*/
static unsigned int v01magic = 0xe3edb0b1; /* v0.1 binary: "cm01" + 0x80808080 */
static unsigned int v01swap  = 0xb1b0ede3; /* v0.1 binary, byteswapped         */

/* Magic tags for cm structure fields
 */
#define CMIO_END_DATA     0 
#define CMIO_NAME         1
#define CMIO_ACC          2
#define CMIO_DESC         3
#define CMIO_ALPHABETTYPE 4
#define CMIO_ALPHABETSIZE 5
#define CMIO_NULL         6
#define CMIO_M            7
#define CMIO_STTYPE       8
#define CMIO_NDIDX        9
#define CMIO_STID         10
#define CMIO_CFIRST       11
#define CMIO_CNUM         12
#define CMIO_PLAST        13
#define CMIO_PNUM         14
#define CMIO_NODES        15
#define CMIO_NODEMAP      16
#define CMIO_NDTYPE       17
#define CMIO_T            18
#define CMIO_E            19
#define CMIO_W            20
#define CMIO_ELSELFSC     21

static void write_ascii_cm(FILE *fp, CM_t *cm);
static int  read_ascii_cm(CMFILE *cmf, CM_t **ret_cm);

static void write_binary_cm(FILE *fp, CM_t *cm);
static int  read_binary_cm(CMFILE *cmf, CM_t **ret_cm);
static void tagged_fwrite(int tag, void *ptr, size_t size, size_t nmemb, FILE *fp);
static int  tagged_fread(int expected_tag, char *s, size_t size, size_t nmemb, FILE *fp);
static void tagged_bin_string_write(int tag, char *s, FILE *fp);
static int  tagged_bin_string_read(int expected_tag, char **ret_s, FILE *fp);

static char *prob2ascii(float p, float null);
static float ascii2prob(char *s, float null);

/* Function:  CMFileOpen()
 * Incept:    SRE, Tue Aug 13 10:25:33 2002 [St. Louis]
 *
 * Purpose:   Opens a CMFILE for reading. This might be a single CM,
 *            or might contain several CMs, or might even be a CM
 *            library with an SSI index; also, it might be in either
 *            binary or ascii format, which we autodetect.
 *            
 *            Looks first in current directory. If cmfile is not
 *            found there, starts looking in the list of directories
 *            in the colon-delimited env string.
 *
 * Args:      cmfile  - name of file to open
 *            env     - NULL, or name of an env variable (e.g. "INFERNALDB"),
 *                      which gives colon-delimited list of directory path(s) 
 *                      to Infernal CM database(s).
 *
 * Returns:   ptr to open CMFILE, or NULL on failure to open.
 *
 * Xref:      STL6 p.108
 */
CMFILE *
CMFileOpen(char *cmfile, char *env)
{
  CMFILE       *cmf;
  unsigned int  magic;
  char         *ssifile;
  char         *dir;
  char          buf[512];

  /* Allocate the CMFILE, and initialize.
   */
  cmf            = MallocOrDie(sizeof(CMFILE));
  cmf->f         = NULL;
  cmf->ssi       = NULL;
  cmf->is_binary = FALSE;
  cmf->byteswap  = FALSE;

  /* Open the file. 
   * Construct ssifile name we'll try to open.
   * Determine recommended index mode (used if we're building an SSI index).
   * Open in mode "r" even if it's binary, not "rb", because we only 
   * guarantee POSIX compatibility, not general ANSI C.
   */
  if ((cmf->f = fopen(cmfile, "r")) != NULL)
    {
      ssifile = MallocOrDie(sizeof(char) * (strlen(cmfile) + 5));
      sprintf(ssifile, "%s.ssi", cmfile);
      
      if ((cmf->mode = SSIRecommendMode(cmfile)) == -1)
	Die("SSIRecommendMode() failed");
    }
  else if ((cmf->f = EnvFileOpen(cmfile, env, &dir)) != NULL)
    {
      char *full;
      full    = FileConcat(dir, cmfile);

      ssifile = MallocOrDie(sizeof(char) * (strlen(full) + strlen(cmfile) + 5));
      sprintf(ssifile, "%s.ssi", full);

      if ((cmf->mode = SSIRecommendMode(full)) == -1)
	Die("SSIRecommendMode() failed unexpectedly");

      free(full);
      free(dir);
    }
  else return NULL;
  
  /* Open the SSI index file, if it exists; if it doesn't,
   * cmf->ssi stays NULL.
   */
  SSIOpen(ssifile, &(cmf->ssi));
  free(ssifile);

  /* Now initialize the disk offset; though it's technically
   * undefined... cmf->offset is the offset of the *last*
   * CM read, so the API only guarantees it's valid after a
   * call to CMFileRead() ... but make it a valid offset 0 
   * anyway. Since the offset is an opaque type, you can't
   * just set it to a number.
   */
  if (SSIGetFilePosition(cmf->f, cmf->mode, &(cmf->offset)) != 0)
    Die("SSIGetFilePosition() failed unexpectedly");
 
  /* Peek at the first 4 bytes to see if it's a binary file.
   */
  if (! fread((char *) &magic, sizeof(unsigned int), 1, cmf->f)) {
    CMFileClose(cmf);
    return NULL;
  }
  rewind(cmf->f);

  /* If the magic number matches one of our binary codes,
   * set the appropriate stuff and return success. Else,
   * fall through to ASCII file tests.
   */
  if (magic == v01magic) {
    cmf->is_binary = TRUE;
    return cmf;
  } else if (magic == v01swap) {
    cmf->is_binary = TRUE;
    cmf->byteswap  = TRUE;
    return cmf;
  } else if (magic & 0x80000000) 
    Die("\
%s appears to be a binary file but the format is not recognized.\n\
It may be from an Infernal version more recent than yours,\n\
or may be a different kind of binary altogether.\n", cmfile);

  /* Check for ASCII format by peeking at first word,
   * and rewind (again!)
   */
  if (fgets(buf, 512, cmf->f) == NULL) {
    CMFileClose(cmf);
    return NULL;
  }
  rewind(cmf->f);

  /* If we recognize the ASCII file tag, return successfully.
   */
  if (strncmp("INFERNAL-1", buf, 10) == 0)
    return cmf;
  
  /* If we haven't recognized the file by now, fail.
   */
  CMFileClose(cmf);
  return NULL;
}

/* Function:  CMFileRead()
 * Incept:    SRE, Tue Aug 13 11:27:55 2002 [St. Louis]
 *
 * Purpose:   Read the next CM in the open file.
 *            Sets the offset in the CMFILE structure to
 *            the offset to the start of this CM.
 *
 * Args:      cmf    - open CMFILE, positioned at start of a CM
 *            ret_cm - RETURN: cm, or NULL on any parsing failure
 *
 * Returns:   1 on success; 0 at EOF.
 *
 * Xref:      STL6 p.108.
 */
int
CMFileRead(CMFILE *cmf, CM_t **ret_cm)
{
  if (SSIGetFilePosition(cmf->f, cmf->mode, &(cmf->offset)) != 0)
    Die("SSIGetFilePosition() failed unexpectedly");
  
  if (cmf->is_binary) return read_binary_cm(cmf, ret_cm);
  else                return read_ascii_cm(cmf, ret_cm);
}

/* Function:  CMFileClose()
 * Incept:    SRE, Tue Aug 13 11:32:58 2002 [St. Louis]
 *
 * Purpose:   Close an open CMFILE.
 *
 * Xref:      STL6 p.108
 */
void
CMFileClose(CMFILE *cmf)
{
  if (cmf->f   != NULL) { fclose(cmf->f);     cmf->f   = NULL; }
  if (cmf->ssi != NULL) { SSIClose(cmf->ssi); cmf->ssi = NULL; }
  free(cmf);
}


/* Function:  CMFileRewind(), CMFilePositionByIndex(), CMFilePositionByKey()
 * Incept:    SRE, Tue Aug 13 11:34:51 2002 [St. Louis]
 *
 * Purpose:   File positioning functions; move to the first CM,
 *            CM #idx (0..ncm-1), or the CM with a given name or accession,
 *            respectively. Return 1 on success, 0 on failure.
 *
 * Xref:      STL6 p.108
 */
void
CMFileRewind(CMFILE *cmf)
{
  rewind(cmf->f);
}
int
CMFilePositionByKey(CMFILE *cmf, char *key)
{
  SSIOFFSET  offset;		/* offset in cmfile, from SSI */
  int        fh;		/* ignored.                   */

  if (cmf->ssi == NULL) return 0;
  if (SSIGetOffsetByName(cmf->ssi, key, &fh, &offset) != 0) return 0;
  if (SSISetFilePosition(cmf->f, &offset) != 0) return 0;
  return 1;
} 
int 
CMFilePositionByIndex(CMFILE *cmf, int idx)
{				/* idx runs from 0..ncm-1 */
  int        fh;		/* file handle is ignored; only one CM file */
  SSIOFFSET  offset;		/* file position of CM */

  if (cmf->ssi == NULL) return 0;
  if (SSIGetOffsetByNumber(cmf->ssi, idx, &fh, &offset) != 0) return 0;
  if (SSISetFilePosition(cmf->f, &offset) != 0) return 0;
  return 1;
}


/* Function:  CMFileWrite()
 * Incept:    SRE, Tue Aug 13 11:41:12 2002 [St. Louis]
 *
 * Purpose:   Write a CM to an open FILE.
 *            If do_binary is TRUE, use binary format; else flatfile.
 * Xref:      STL6 p.108.
 */
void
CMFileWrite(FILE *fp, CM_t *cm, int do_binary)
{
  if (do_binary) write_binary_cm(fp, cm);
  else           write_ascii_cm(fp, cm);
}
		   

/* Function:  write_ascii_cm()
 * Incept:    SRE, Tue Aug 13 11:45:43 2002 [St. Louis]
 *
 * Purpose:   Write a CM in flatfile format.
 * Xref:      STL6 p.108
 */
static void
write_ascii_cm(FILE *fp, CM_t *cm)
{
  int v,x,y,nd;

  fprintf(fp, "INFERNAL-1 [%s]\n", PACKAGE_VERSION);

  fprintf(fp, "NAME   %s\n", cm->name);
  if (cm->acc  != NULL)  fprintf(fp, "ACC    %s\n", cm->acc);
  if (cm->desc != NULL)  fprintf(fp, "DESC   %s\n", cm->desc);
  fprintf(fp, "STATES %d\n", cm->M);
  fprintf(fp, "NODES  %d\n", cm->nodes);
  /* EPN 08.18.05 */
  fprintf(fp, "W      %d\n", cm->W);
  /* EPN 11.15.05 */
  fprintf(fp, "el_selfsc %f\n", cm->el_selfsc);

  fputs("NULL  ", fp);
  for (x = 0; x < Alphabet_size; x++)
    fprintf(fp, "%6s ", prob2ascii(cm->null[x], 1/(float)(Alphabet_size)));
  fputs("\n", fp);

  fputs("MODEL:\n", fp);
  for (v = 0; v < cm->M; v++) 
    {
      nd = cm->ndidx[v];

      /* Node line.
       */
      if (cm->nodemap[nd] == v) 
	fprintf(fp, "\t\t\t\t[ %-4s %4d ]\n", Nodetype(cm->ndtype[nd]), nd);

      /* State line, w/ parents, children, and transitions
       */
      fprintf(fp, "    %2s %5d %5d %1d %5d %5d ", 
	      Statetype(cm->sttype[v]), v, 
	      cm->plast[v], cm->pnum[v],
	      cm->cfirst[v], cm->cnum[v]);
      if (cm->sttype[v] != B_st)
	for (x = 0; x < cm->cnum[v]; x++)
	  fprintf(fp, "%7s ", prob2ascii(cm->t[v][x], 1.));
      else x = 0;
      for (; x < 6; x++)
	fprintf(fp, "%7s ", "");
      
      /* Emission line
       */
      if (cm->sttype[v] == MP_st) 
	{
	  for (x = 0; x < Alphabet_size; x++)
	    for (y = 0; y < Alphabet_size; y++)
	      fprintf(fp, "%6s ", prob2ascii(cm->e[v][x*Alphabet_size+y], cm->null[x]*cm->null[y]));
	}
      else if (cm->sttype[v] == ML_st || cm->sttype[v] == MR_st || cm->sttype[v] == IL_st || cm->sttype[v] == IR_st)
	{
	  for (x = 0; x < Alphabet_size; x++)
	    fprintf(fp, "%6s ", prob2ascii(cm->e[v][x], cm->null[x]));
	}
      fputs("\n", fp);
    }
  fputs("//\n", fp);
} 

static int  
read_ascii_cm(CMFILE *cmf, CM_t **ret_cm)
{
  CM_t   *cm;
  char   *buf;
  int     n;			/* length of buf */
  char   *s;
  int     M,N;			/* number of states, nodes in model */
  int     v,x,y,nd;		/* counters for states, events, nodes */
  char   *tok;
  int     toklen;

  cm  = NULL;
  buf = NULL;
  n   = 0;
  if (feof(cmf->f) || sre_fgets(&buf, &n, cmf->f) == NULL) return 0;
  if (strncmp(buf, "INFERNAL-1", 10) != 0)                 goto FAILURE;

  /* Parse the header information
   * These are all tag/value. 
   * Ignore unknown tags (forward compatibility). 
   */
  cm = CreateCMShell();
  M  = N = -1;
  while (sre_fgets(&buf, &n, cmf->f) != NULL) 
    {
      s   = buf;
      if ((tok = sre_strtok(&s, " \t\n", &toklen)) == NULL) goto FAILURE;
      if (strcmp(tok, "NAME") == 0) 
	{
	  if ((tok = sre_strtok(&s, " \t\n", &toklen)) == NULL) goto FAILURE;
	  cm->name = sre_strdup(tok, toklen);
	}
      else if (strcmp(tok, "ACC") == 0) 
	{
	  if ((tok = sre_strtok(&s, " \t\n", &toklen)) == NULL) goto FAILURE;
	  cm->acc = sre_strdup(tok, toklen);
	}
      else if (strcmp(tok, "DESC") == 0) 
	{
	  if ((tok = sre_strtok(&s, " \t\n", &toklen)) == NULL) goto FAILURE;
	  cm->desc = sre_strdup(tok, toklen);
	}
      else if (strcmp(tok, "STATES") == 0) 
	{
	  if ((tok = sre_strtok(&s, " \t\n", &toklen)) == NULL) goto FAILURE;
	  M = atoi(tok);
	}
      else if (strcmp(tok, "NODES") == 0) 
	{
	  if ((tok = sre_strtok(&s, " \t\n", &toklen)) == NULL) goto FAILURE;
	  N = atoi(tok);
	}
      else if (strcmp(tok, "NULL") == 0) 
	{
	  for (x = 0; x < Alphabet_size; x++)
	    {
	      if ((tok = sre_strtok(&s, " \t\n", &toklen)) == NULL) goto FAILURE;
	      cm->null[x] = ascii2prob(tok, (1./(float)Alphabet_size));
	    }
	}
      /* EPN 08.18.05 */
      else if (strcmp(tok, "W") == 0) 
	{
	  if ((tok = sre_strtok(&s, " \t\n", &toklen)) == NULL) goto FAILURE;
	  cm->W = atoi(tok);
	}
      /* EPN 11.15.05 */
      else if (strcmp(tok, "el_selfsc") == 0) 
	{
	  if ((tok = sre_strtok(&s, " \t\n", &toklen)) == NULL) goto FAILURE;
	  cm->el_selfsc = atof(tok);
	}
      else if (strcmp(tok, "MODEL:") == 0)
	break;
    }

  /* Done reading the header information.
   * Check that everything is ok and mandatory info is present before moving on.
   */
  if (feof(cmf->f))       goto FAILURE;
  if (M < 1)              goto FAILURE;
  if (N < 1)              goto FAILURE;
  if (cm->name == NULL)   goto FAILURE;

  /* Main model section. 
   */
  CreateCMBody(cm, N, M);
  CMZero(cm);
  nd = -1;
  for (v = 0; v < cm->M; v++)
    {
      if (sre_fgets(&buf, &n, cmf->f) == NULL) goto FAILURE;
      s = buf;
      if ((tok = sre_strtok(&s, " \t\n", &toklen)) == NULL) goto FAILURE;      
      
      /* Ah, a node line. Process it and get the following line.
       */
      if (*tok == '[') 
	{
	  if ((tok = sre_strtok(&s, " \t\n", &toklen)) == NULL) goto FAILURE;      
	  if ((x = NodeCode(tok)) == -1)                        goto FAILURE;
	  if ((tok = sre_strtok(&s, " \t\n", &toklen)) == NULL) goto FAILURE;      
	  if (!IsInt(tok))                                      goto FAILURE;
	  nd = atoi(tok);
	  cm->ndtype[nd]  = x;
	  cm->nodemap[nd] = v;

	  if (sre_fgets(&buf, &n, cmf->f) == NULL)              goto FAILURE;
	  s = buf;
	  if ((tok = sre_strtok(&s, " \t\n", &toklen)) == NULL) goto FAILURE;      
	}

      /* Process state line.
       */
      if ((cm->sttype[v] = StateCode(tok)) == -1)           goto FAILURE;
      if ((tok = sre_strtok(&s, " \t\n", &toklen)) == NULL) goto FAILURE;      
      if (! IsInt(tok))                                     goto FAILURE;
      if (atoi(tok) != v)                                   goto FAILURE;
      if ((tok = sre_strtok(&s, " \t\n", &toklen)) == NULL) goto FAILURE;      
      if (! IsInt(tok))                                     goto FAILURE;
      cm->plast[v] = atoi(tok);
      if ((tok = sre_strtok(&s, " \t\n", &toklen)) == NULL) goto FAILURE;      
      if (! IsInt(tok))                                     goto FAILURE;
      cm->pnum[v] = atoi(tok);
      if ((tok = sre_strtok(&s, " \t\n", &toklen)) == NULL) goto FAILURE;      
      if (! IsInt(tok))                                     goto FAILURE;
      cm->cfirst[v] = atoi(tok);
      if ((tok = sre_strtok(&s, " \t\n", &toklen)) == NULL) goto FAILURE;      
      if (! IsInt(tok))                                     goto FAILURE;
      cm->cnum[v] = atoi(tok);
				/* Transition probabilities. */
      if (cm->sttype[v] != B_st) 
	{
	  for (x = 0; x < cm->cnum[v]; x++)
	    {
	      if ((tok = sre_strtok(&s, " \t\n", &toklen)) == NULL) goto FAILURE;      
	      if (! IsReal(tok) && *tok != '*')                      goto FAILURE;
	      cm->t[v][x] = ascii2prob(tok, 1.);
	    }
	}
				/* Emission probabilities. */
      if (cm->sttype[v] == ML_st || cm->sttype[v] == MR_st ||
	  cm->sttype[v] == IL_st || cm->sttype[v] == IR_st)
	{
	  for (x = 0; x < Alphabet_size; x++)
	    {
	      if ((tok = sre_strtok(&s, " \t\n", &toklen)) == NULL) goto FAILURE;      
	      if (! IsReal(tok) && *tok != '*')                     goto FAILURE;
	      cm->e[v][x] = ascii2prob(tok, cm->null[x]);
	    }
	}
      else if (cm->sttype[v] == MP_st) 
	{
	  for (x = 0; x < Alphabet_size; x++)
	    for (y = 0; y < Alphabet_size; y++)
	      {
		if ((tok = sre_strtok(&s, " \t\n", &toklen)) == NULL) goto FAILURE;      
		if (! IsReal(tok) && *tok != '*')                     goto FAILURE;
		cm->e[v][x*Alphabet_size+y] = ascii2prob(tok, cm->null[x]*cm->null[y]);
	      }
	} 

      cm->ndidx[v] = nd;
      cm->stid[v]  = DeriveUniqueStateCode(cm->ndtype[nd], cm->sttype[v]);
    } /* end of loop over states */

  /* Advance to record separator
   */
  while (sre_fgets(&buf, &n, cmf->f) != NULL) 
    if (strncmp(buf, "//", 2) == 0) 
      break;

  /* EPN 10.29.06 Remove the sole source of CM ambiguities. Find and detach insert states
   *              that are 1 state before an END_E.  */
  cm_find_and_detach_dual_inserts(cm, 
				  FALSE, /* Don't check END_E-1 states have 0 counts, they may not if 
					  * an old version (0.7 or earlier) of cmbuild was used, or  
					  * cmbuild --nodetach  was used to build the CM  */
				  TRUE); /* Detach the states by setting trans probs into them as 0.0   */

  /* Success.
   * Renormalize the CM, build the CP9 HMM and return.
   */
  CMRenormalize(cm);

  /* Build the CM Plan 9 HMM, the HMM and cp9map are allocated in build_cp9_hm() */
  if(!build_cp9_hmm(cm, &(cm->cp9), &(cm->cp9map), FALSE, 0.001, 0))
    Die("Couldn't build a CP9 HMM from the CM\n");



  if (buf != NULL) free(buf);
  *ret_cm = cm;
  return 1;

 FAILURE:
  if (cm != NULL)  FreeCM(cm);
  if (buf != NULL) free(buf);
  *ret_cm = NULL;
  return 1;
}

/* Function: write_binary_cm()
 * Date:     SRE, Thu Aug  3 12:05:30 2000 [St. Louis]
 *
 * Purpose:  Write a CM in binary format.
 *
 */
static void
write_binary_cm(FILE *fp, CM_t *cm)
{
  int v;

  fwrite((char *) &(v01magic), sizeof(unsigned int), 1, fp);

  /* These have to go first, so we know how big of a CM
   * to allocate when we go to read a file.
   */
  tagged_fwrite(CMIO_M,            &cm->M,          sizeof(int),   1, fp);
  tagged_fwrite(CMIO_NODES,        &cm->nodes,      sizeof(int),   1, fp);  

  tagged_bin_string_write(CMIO_NAME, cm->name,  fp);
  tagged_bin_string_write(CMIO_ACC,  cm->acc,   fp);
  tagged_bin_string_write(CMIO_DESC, cm->desc,  fp);

  tagged_fwrite(CMIO_ALPHABETTYPE, &Alphabet_type, sizeof(int),   1, fp);
  tagged_fwrite(CMIO_ALPHABETSIZE, &Alphabet_size, sizeof(int),   1, fp);
  tagged_fwrite(CMIO_NULL,         cm->null,       sizeof(float), Alphabet_size, fp);
  tagged_fwrite(CMIO_STTYPE,       cm->sttype,     sizeof(char),  cm->M, fp);
  tagged_fwrite(CMIO_NDIDX,        cm->ndidx,      sizeof(int),   cm->M, fp);  
  tagged_fwrite(CMIO_STID,         cm->stid,       sizeof(char),  cm->M, fp);  
  tagged_fwrite(CMIO_CFIRST,       cm->cfirst,     sizeof(int),   cm->M, fp); 
  tagged_fwrite(CMIO_CNUM,         cm->cnum,       sizeof(int),   cm->M, fp);  
  tagged_fwrite(CMIO_PLAST,        cm->plast,      sizeof(int),   cm->M, fp); 
  tagged_fwrite(CMIO_PNUM,         cm->pnum,       sizeof(int),   cm->M, fp);  
  tagged_fwrite(CMIO_NODEMAP,      cm->nodemap,    sizeof(int),   cm->nodes, fp);
  tagged_fwrite(CMIO_NDTYPE,       cm->ndtype,     sizeof(char),  cm->nodes, fp);
  /* EPN 08.18.05 */
  tagged_fwrite(CMIO_W,           &cm->W,          sizeof(int),    1, fp);  
  /* EPN 11.15.05 */
  tagged_fwrite(CMIO_ELSELFSC,    &cm->el_selfsc,   sizeof(float),  1, fp);  

  for (v = 0; v < cm->M; v++) {
    tagged_fwrite(CMIO_T, cm->t[v], sizeof(float), MAXCONNECT, fp);
    tagged_fwrite(CMIO_E, cm->e[v], sizeof(float), Alphabet_size*Alphabet_size, fp);
  }

  tagged_fwrite(CMIO_END_DATA, NULL, 0, 0, fp);

  /* Note: begin, end, and flags not written out. Local alignment is
   * run-time configuration right now.
   */
}


/* Function: read_binary_cm()
 * Date:     SRE, Thu Aug  3 13:39:09 2000 [St. Louis]
 *
 * Purpose:  Read a CM from disk.
 */
static int
read_binary_cm(CMFILE *cmf, CM_t **ret_cm)
{
  FILE         *fp;
  CM_t         *cm;
  unsigned int  magic;
  int           M;
  int           nodes;
  int           atype;
  int           asize;
  int           v;

  cm = NULL;
  fp = cmf->f;
  if (feof(fp)) return 0;
  if (! fread((char *) &magic, sizeof(unsigned int), 1, fp)) return 0;
  if (magic != v01magic) goto FAILURE;
  
  if (! tagged_fread(CMIO_M,     (void *) &M,     sizeof(int), 1, fp)) goto FAILURE;
  if (! tagged_fread(CMIO_NODES, (void *) &nodes, sizeof(int), 1, fp)) goto FAILURE;
  cm = CreateCM(nodes, M);

  if (! tagged_bin_string_read(CMIO_NAME, &(cm->name),  fp)) goto FAILURE;
  if (! tagged_bin_string_read(CMIO_ACC,  &(cm->acc),   fp)) goto FAILURE;
  if (! tagged_bin_string_read(CMIO_DESC, &(cm->desc),  fp)) goto FAILURE;
  
  if (! tagged_fread(CMIO_ALPHABETTYPE, (void *) &atype,         sizeof(int),   1, fp))             goto FAILURE;
  if (! tagged_fread(CMIO_ALPHABETSIZE, (void *) &asize,         sizeof(int),   1, fp))             goto FAILURE;
  if (! tagged_fread(CMIO_NULL,         (void *) cm->null,       sizeof(float), Alphabet_size, fp)) goto FAILURE;
  if (! tagged_fread(CMIO_STTYPE,       (void *) cm->sttype,     sizeof(char),  cm->M, fp))         goto FAILURE;
  if (! tagged_fread(CMIO_NDIDX,        (void *) cm->ndidx,      sizeof(int),   cm->M, fp))         goto FAILURE;  
  if (! tagged_fread(CMIO_STID,         (void *) cm->stid,       sizeof(char),  cm->M, fp))         goto FAILURE;  
  if (! tagged_fread(CMIO_CFIRST,       (void *) cm->cfirst,     sizeof(int),   cm->M, fp))         goto FAILURE; 
  if (! tagged_fread(CMIO_CNUM,         (void *) cm->cnum,       sizeof(int),   cm->M, fp))         goto FAILURE;
  if (! tagged_fread(CMIO_PLAST,        (void *) cm->plast,      sizeof(int),   cm->M, fp))         goto FAILURE; 
  if (! tagged_fread(CMIO_PNUM,         (void *) cm->pnum,       sizeof(int),   cm->M, fp))         goto FAILURE;  
  if (! tagged_fread(CMIO_NODEMAP,      (void *) cm->nodemap,    sizeof(int),   cm->nodes, fp))     goto FAILURE;
  if (! tagged_fread(CMIO_NDTYPE,       (void *) cm->ndtype,     sizeof(char),  cm->nodes, fp))     goto FAILURE;
  /* EPN 08.18.05 */
  if (! tagged_fread(CMIO_W,            (void *) &(cm->W),       sizeof(int),   1,         fp))     goto FAILURE;

  if (! tagged_fread(CMIO_ELSELFSC,     (void *) &(cm->el_selfsc), sizeof(float),1,        fp))     goto FAILURE;
  
  for (v = 0; v < cm->M; v++) {
    if (! tagged_fread(CMIO_T, (void *) cm->t[v], sizeof(float), MAXCONNECT, fp)) goto FAILURE;
    if (! tagged_fread(CMIO_E, (void *) cm->e[v], sizeof(float), Alphabet_size*Alphabet_size, fp)) goto FAILURE;
  }
  if (! tagged_fread(CMIO_END_DATA, (void *) NULL, 0, 0, fp)) goto FAILURE;

  /* EPN 10.29.06 Remove the sole source of CM ambiguities. Find and detach insert states
   *              that are 1 state before an END_E.  */
  cm_find_and_detach_dual_inserts(cm, 
				  FALSE, /* Don't check END_E-1 states have 0 counts, they may not if 
					  * an old version (0.7 or earlier) of cmbuild was used, or  
					  * cmbuild --nodetach  was used to build the CM  */
				  TRUE); /* Detach the states by setting trans probs into them as 0.0   */

  /* Build the CM Plan 9 HMM, the HMM and cp9map are allocated in build_cp9_hm() */
  if(!build_cp9_hmm(cm, &(cm->cp9), &(cm->cp9map), FALSE, 0.001, 0))
    Die("Couldn't build a CP9 HMM from the CM\n");

  /* EPN 10.29.06 Noticed there's no CMRenormalize() call here (for speed?), didn't add one 
     CMRenormalize(cm);
   */
  *ret_cm = cm;
  return 1;

 FAILURE:
  if (cm != NULL) FreeCM(cm);
  *ret_cm = NULL;
  return 1;
}

static void
tagged_fwrite(int tag, void *ptr, size_t size, size_t nmemb, FILE *fp)
{
  int nbytes;

  nbytes = size * nmemb;
  fwrite(&tag,    sizeof(int), 1,     fp);
  fwrite(&nbytes, sizeof(int), 1,     fp);
  if (nbytes != 0 && ptr != NULL) 
    fwrite(ptr,   size,        nmemb, fp);
}
static int
tagged_fread(int expected_tag, char *s, size_t size, size_t nmemb, FILE *fp)
{
  int tag;
  int nbytes;
  
  fread(&tag,    sizeof(int), 1,     fp); if (tag    != expected_tag) return 0;
  fread(&nbytes, sizeof(int), 1,     fp); if (nbytes != (size*nmemb)) return 0;
  if (nbytes != 0)
    fread(s, size, nmemb, fp); 
  return 1;
}
static void
tagged_bin_string_write(int tag, char *s, FILE *fp)
{
  int len;
  if (s != NULL) {
    len = strlen(s);
    tagged_fwrite(tag, s, sizeof(char), len, fp);
  } else {
    tagged_fwrite(tag, NULL, 0, 0, fp);
  }
}
static int
tagged_bin_string_read(int expected_tag, char **ret_s, FILE *fp)
{
  int tag;
  int nbytes;
  char *s;

  fread(&tag, sizeof(int), 1, fp);
  if (tag != expected_tag) return 0;
  fread(&nbytes, sizeof(int), 1, fp);
  if (nbytes > 0) {
    s = MallocOrDie(sizeof(char) * (nbytes+1));
    s[nbytes] = '\0';
    fread(s, sizeof(char), nbytes, fp);
  } else s = NULL; 
  *ret_s = s;
  return 1;
}
/*****************************************************************
 * Some miscellaneous utility functions
 *****************************************************************/

/* Function: prob2ascii()
 * 
 * Purpose:  Format a probability for output to an ASCII save
 *           file. Returns a ptr to a static internal buffer.
 *              
 */
static char *
prob2ascii(float p, float null)
{
  static char buffer[32];

  if (p == 0.0) return "*";
  sprintf(buffer, "%.3f", sreLOG2(p/null));
  return buffer;
}


/* Function: ascii2prob()
 * 
 * Purpose:  Convert a saved string back to a probability.
 */
static float
ascii2prob(char *s, float null)
{
  return (*s == '*') ? 0. : exp(atof(s)/1.44269504)*null;
}

