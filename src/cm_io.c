/* cm_io.c
 * SRE, Thu Aug  3 11:53:34 2000 [St. Louis]
 * SVN $Id$
 * 
 * Input/output of covariance models.
 * 
 *****************************************************************
 * @LICENSE@
 ***************************************************************** 
 */ 

#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_ssi.h"

#include "funcs.h"
#include "structs.h"

static int is_integer(char *s);
static int is_real(char *s);

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
#define CMIO_NPART        22
#define CMIO_PARTS        23
#define CMIO_PARTE        24
#define CMIO_EXPLAMBDA    25
#define CMIO_EXPMUE       26
#define CMIO_EXPMUO       27
#define CMIO_EXPDBSIZE    28
#define CMIO_EXPNHITS     29
#define CMIO_EXPTAILP     30
#define CMIO_FTHRNCUT     31
#define CMIO_FTHRF        32
#define CMIO_FTHRN        33
#define CMIO_FTHRDB       34
#define CMIO_FTHRBETA     35
#define CMIO_FTHRUSEQDB   36
#define CMIO_FTHRABTS     37
#define CMIO_FTHRCMECUT   38
#define CMIO_FTHRFWDECUT  39
#define CMIO_HASEXP       40
#define CMIO_HASFILTER    41
#define CMIO_ABCTYPE      42
#define CMIO_HASGA        43
#define CMIO_HASTC        44
#define CMIO_HASNC        45
#define CMIO_GA           46
#define CMIO_TC           47
#define CMIO_NC           48
#define CMIO_BCOM         49
#define CMIO_BDATE        50
#define CMIO_CCOM         51
#define CMIO_CDATE        52
#define CMIO_NSEQ         53
#define CMIO_EFFNSEQ      54
#define CMIO_CLEN         55
#define CMIO_WBETA        56

static int  write_ascii_cm(FILE *fp, CM_t *cm, char *errbuf);
static int  read_ascii_cm(CMFILE *cmf, char *errbuf, ESL_ALPHABET **ret_abc, CM_t **ret_cm);
static int  write_binary_cm(FILE *fp, CM_t *cm, char *errbuf);
static int  read_binary_cm(CMFILE *cmf, char *errbuf, ESL_ALPHABET **ret_abc, CM_t **ret_cm);
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
  int           status;
  CMFILE       *cmf;
  unsigned int  magic;
  char         *ssifile = NULL;	/* constructed name of SSI index file             */
  char         *envfile = NULL;	/* full path to filename after using environment  */
  char          buf[512];
  int           n = strlen(cmfile);

  /* Allocate the CMFILE, and initialize.
   */
  ESL_ALLOC(cmf, sizeof(CMFILE));
  cmf->f         = NULL;
  cmf->fname     = NULL;
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
      if ((status = esl_FileNewSuffix(cmfile, "ssi", &ssifile)) != eslOK) goto ERROR;
      if ((status = esl_strdup(cmfile, n, &(cmf->fname)))       != eslOK) goto ERROR;
    }
  else if (esl_FileEnvOpen(cmfile, env, &(cmf->f), &envfile) == eslOK)
    {
      if ((status = esl_FileNewSuffix(envfile, "ssi", &ssifile)) != eslOK) goto ERROR;
      if ((status = esl_strdup(envfile, -1, &(cmf->fname)))      != eslOK) goto ERROR;
    }
  else
    { status = eslENOTFOUND; goto ERROR; }

  /* Attempt to open the ssi index file. cmf->ssi silently stays NULL if the ssifile isn't found. */
  if (ssifile != NULL) esl_ssi_Open(ssifile, &(cmf->ssi));
  if (envfile != NULL) free(envfile);
  if (ssifile != NULL) free(ssifile);

  /* Now initialize the disk offset; though it's technically
   * undefined... cmf->offset is the offset of the *last*
   * CM read, so the API only guarantees it's valid after a
   * call to CMFileRead() ... but make it a valid offset 0 
   * anyway. Since the offset is an opaque type, you can't
   * just set it to a number.
   */
  cmf->offset = ftello(cmf->f); 

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
    cm_Fail("\
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

 ERROR:
  return NULL; 
}

/* Function:  CMFileRead()
 * Incept:    SRE, Tue Aug 13 11:27:55 2002 [St. Louis]
 *
 * Purpose:   Read the next CM in the open file.
 *            Sets the offset in the CMFILE structure to
 *            the offset to the start of this CM.
 *
 *            From HMMER3:
 *            Caller may or may not already know what alphabet the HMM
 *            is expected to be in.  A reference to the pointer to the
 *            current alphabet is passed in <*ret_abc>. If the alphabet
 *            is unknown, pass <*ret_abc = NULL>, and when the
 *            new CM is read, an appropriate new alphabet object is
 *            allocated and passed back to the caller in <*ret_abc>.
 *            If the alphabet is already known, <ret_abc> points to
 *            that object ptr, and the new HMM's alphabet type is
 *            verified to agree with it. This mechanism allows an
 *            application to let the first HMM determine the alphabet
 *            type for the application, while still keeping the
 *            alphabet under the application's scope of control.
 *
 * Args:      cmf    - open CMFILE, positioned at start of a CM
 *            ret_abc- alphabet (see above)
 *            ret_cm - RETURN: cm, or NULL on any parsing failure
 *
 * Returns:   eslOK on success; eslEOF at EOF.
 * Throws:    eslEFORMAT if format of file is incorrect (pre 1.0 cmfile for example)
 *            errbuf is filled with error message.
 *
 * Xref:      STL6 p.108.
 */
int
CMFileRead(CMFILE *cmf, char *errbuf, ESL_ALPHABET **ret_abc, CM_t **ret_cm)
{
  cmf->offset = ftello(cmf->f); 
  if (cmf->is_binary) return read_binary_cm(cmf, errbuf, ret_abc, ret_cm);
  else                return read_ascii_cm(cmf, errbuf, ret_abc, ret_cm);
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
  if (cmf->f     != NULL) { fclose(cmf->f);     cmf->f   = NULL; }
  if (cmf->fname != NULL)   free(cmf->fname); 
  if (cmf->ssi   != NULL) { esl_ssi_Close(cmf->ssi); cmf->ssi = NULL; }
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
  uint16_t fh;
  off_t    offset;
  int      status;

  if (cmf->ssi == NULL) ESL_EXCEPTION(eslEINVAL, "CMFilePositionByKey(): cmf->ssi is NULL");
  if ((status = esl_ssi_FindName(cmf->ssi, key, &fh, &offset, NULL, NULL)) != eslOK) return status;
  if (fseeko(cmf->f, offset, SEEK_SET) != 0)    ESL_EXCEPTION(eslESYS, "fseek failed");
  return eslOK;
} 
int 
CMFilePositionByIndex(CMFILE *cmf, int64_t idx)
{				/* idx runs from 0..ncm-1 */
  uint16_t fh;
  off_t    offset;
  int      status;

  if (cmf->ssi == NULL) ESL_EXCEPTION(eslEINVAL, "CMFilePositionByIndex(): cmf->ssi is NULL");
  if ((status = esl_ssi_FindNumber(cmf->ssi, idx, &fh, &offset, NULL, NULL, NULL)) != eslOK) return status;
  if (fseeko(cmf->f, offset, SEEK_SET) != 0)    ESL_EXCEPTION(eslESYS, "fseek failed");
  return eslOK;
}


/* Function:  CMFileWrite()
 * Incept:    SRE, Tue Aug 13 11:41:12 2002 [St. Louis]
 *
 * Purpose:   Write a CM to an open FILE.
 *            If do_binary is TRUE, use binary format; else flatfile.
 * Xref:      STL6 p.108.
 *
 * Returns: eslOK on success.
 *          eslEINCOMPAT on contract violation (if a mandatory part of the CM is invalid)
 */
int 
CMFileWrite(FILE *fp, CM_t *cm, int do_binary, char *errbuf)
{
  /* contract checks */
  if (cm->name          == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "CMFileWrite(), cm->name is NULL.");
  if (cm->comlog        == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "CMFileWrite(), cm->comlog is NULL.");
  if (cm->comlog->bcom  == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "CMFileWrite(), cm->comlog->bcom is NULL.");
  if (cm->comlog->bdate == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "CMFileWrite(), cm->comlog->bdate is NULL.");
  if((cm->flags & CMH_LOCAL_BEGIN) && (cm->flags & CMH_LOCAL_END)) ESL_FAIL(eslEINCOMPAT, errbuf, "CMFileWrite(), CMH_LOCAL_BEGIN and CMH_LOCAL_END flags are up.");
  if (cm->flags & CMH_LOCAL_BEGIN) ESL_FAIL(eslEINCOMPAT, errbuf, "CMFileWrite(), CMH_LOCAL_BEGIN flag is up.");
  if (cm->flags & CMH_LOCAL_END)   ESL_FAIL(eslEINCOMPAT, errbuf, "CMFileWrite(), CMH_LOCAL_END flag is up.");
  if (do_binary) return write_binary_cm(fp, cm, errbuf);
  else           return write_ascii_cm(fp, cm, errbuf);
}
		   

/* Function:  write_ascii_cm()
 * Incept:    SRE, Tue Aug 13 11:45:43 2002 [St. Louis]
 *
 * Purpose:   Write a CM in flatfile format.
 * Xref:      STL6 p.108
 *
 * Returns: eslOK on success;
 */
static int
write_ascii_cm(FILE *fp, CM_t *cm, char *errbuf)
{
  int v,x,y,nd,i;
  
  fprintf(fp, "INFERNAL-1 [%s]\n", PACKAGE_VERSION);

  fprintf(fp,                          "NAME     %s\n", cm->name);
  if (cm->acc  != NULL)    fprintf(fp, "ACC      %s\n", cm->acc);
  if (cm->desc != NULL)    fprintf(fp, "DESC     %s\n", cm->desc);
  /* Rfam cutoffs */
  if (cm->flags & CMH_GA)  fprintf(fp, "GA       %.2f\n", cm->ga);
  if (cm->flags & CMH_TC)  fprintf(fp, "TC       %.2f\n", cm->tc);
  if (cm->flags & CMH_NC)  fprintf(fp, "NC       %.2f\n", cm->nc);
  fprintf(fp, "STATES   %d\n",   cm->M);
  fprintf(fp, "NODES    %d\n",   cm->nodes);
  fprintf(fp, "ALPHABET %d\n",   cm->abc->type);
  fprintf(fp, "ELSELF   %.8f\n", cm->el_selfsc);
  fprintf(fp, "WBETA    %g\n",   cm->beta_W);
  fprintf(fp, "NSEQ     %d\n",   cm->nseq);
  fprintf(fp, "EFFNSEQ  %.3f\n", cm->eff_nseq);
  fprintf(fp, "CLEN     %d\n",   cm->clen);
  fprintf(fp, "BCOM     %s\n",   cm->comlog->bcom);
  fprintf(fp, "BDATE    %s\n",   cm->comlog->bdate);
  if(cm->comlog->ccom != NULL) fprintf(fp, "CCOM     %s\n", cm->comlog->ccom);
  if(cm->comlog->cdate!= NULL) fprintf(fp, "CDATE    %s\n", cm->comlog->cdate);
  fputs(      "NULL    ", fp);
  for (x = 0; x < cm->abc->K; x++)
    fprintf(fp, "%6s ", prob2ascii(cm->null[x], 1/(float)(cm->abc->K)));
  fputs("\n", fp);

  /* E-value statistics
   */
  int p;
  if (cm->flags & CMH_EXPTAIL_STATS)
    {
      fprintf(fp, "PART     %-3d  ", cm->stats->np);
      for(p = 0; p < cm->stats->np; p++)
	fprintf(fp, "%5d  %5d  ", cm->stats->ps[p], cm->stats->pe[p]);
      fprintf(fp, "\n");
      for(p = 0; p < cm->stats->np; p++)
	{
	  fprintf(fp, "E-LC     %-2d  %10.5f  %10.5f  %10.5f  %10ld  %10d  %.6f\n", 
		  p, cm->stats->expAA[EXP_CM_LC][p]->lambda, cm->stats->expAA[EXP_CM_LC][p]->mu_extrap, cm->stats->expAA[EXP_CM_LC][p]->mu_orig, 
		  cm->stats->expAA[EXP_CM_LC][p]->dbsize, cm->stats->expAA[EXP_CM_LC][p]->nrandhits, cm->stats->expAA[EXP_CM_LC][p]->tailp);
	  fprintf(fp, "E-GC     %-2d  %10.5f  %10.5f  %10.5f  %10ld  %10d  %.6f\n", 
		  p, cm->stats->expAA[EXP_CM_GC][p]->lambda, cm->stats->expAA[EXP_CM_GC][p]->mu_extrap, cm->stats->expAA[EXP_CM_GC][p]->mu_orig, 
		  cm->stats->expAA[EXP_CM_GC][p]->dbsize, cm->stats->expAA[EXP_CM_GC][p]->nrandhits, cm->stats->expAA[EXP_CM_GC][p]->tailp);
	  fprintf(fp, "E-LI     %-2d  %10.5f  %10.5f  %10.5f  %10ld  %10d  %.6f\n", 
		  p, cm->stats->expAA[EXP_CM_LI][p]->lambda, cm->stats->expAA[EXP_CM_LI][p]->mu_extrap, cm->stats->expAA[EXP_CM_LI][p]->mu_orig, 
		  cm->stats->expAA[EXP_CM_LI][p]->dbsize, cm->stats->expAA[EXP_CM_LI][p]->nrandhits, cm->stats->expAA[EXP_CM_LI][p]->tailp);
	  fprintf(fp, "E-GI     %-2d  %10.5f  %10.5f  %10.5f  %10ld  %10d  %.6f\n", 
		  p, cm->stats->expAA[EXP_CM_GI][p]->lambda, cm->stats->expAA[EXP_CM_GI][p]->mu_extrap, cm->stats->expAA[EXP_CM_GI][p]->mu_orig, 
		  cm->stats->expAA[EXP_CM_GI][p]->dbsize, cm->stats->expAA[EXP_CM_GI][p]->nrandhits, cm->stats->expAA[EXP_CM_GI][p]->tailp);
	  fprintf(fp, "E-LV     %-2d  %10.5f  %10.5f  %10.5f  %10ld  %10d  %.6f\n", 
		  p, cm->stats->expAA[EXP_CP9_LV][p]->lambda, cm->stats->expAA[EXP_CP9_LV][p]->mu_extrap, cm->stats->expAA[EXP_CP9_LV][p]->mu_orig, 
		  cm->stats->expAA[EXP_CP9_LV][p]->dbsize, cm->stats->expAA[EXP_CP9_LV][p]->nrandhits, cm->stats->expAA[EXP_CP9_LV][p]->tailp);
	  fprintf(fp, "E-GV     %-2d  %10.5f  %10.5f  %10.5f  %10ld  %10d  %.6f\n", 
		  p, cm->stats->expAA[EXP_CP9_GV][p]->lambda, cm->stats->expAA[EXP_CP9_GV][p]->mu_extrap, cm->stats->expAA[EXP_CP9_GV][p]->mu_orig, 
		  cm->stats->expAA[EXP_CP9_GV][p]->dbsize, cm->stats->expAA[EXP_CP9_GV][p]->nrandhits, cm->stats->expAA[EXP_CP9_GV][p]->tailp);
	  fprintf(fp, "E-LF     %-2d  %10.5f  %10.5f  %10.5f  %10ld  %10d  %.6f\n", 
		  p, cm->stats->expAA[EXP_CP9_LF][p]->lambda, cm->stats->expAA[EXP_CP9_LF][p]->mu_extrap, cm->stats->expAA[EXP_CP9_LF][p]->mu_orig, 
		  cm->stats->expAA[EXP_CP9_LF][p]->dbsize, cm->stats->expAA[EXP_CP9_LF][p]->nrandhits, cm->stats->expAA[EXP_CP9_LF][p]->tailp);
	  fprintf(fp, "E-GF     %-2d  %10.5f  %10.5f  %10.5f  %10ld  %10d  %.6f\n", 
		  p, cm->stats->expAA[EXP_CP9_GF][p]->lambda, cm->stats->expAA[EXP_CP9_GF][p]->mu_extrap, cm->stats->expAA[EXP_CP9_GF][p]->mu_orig, 
		  cm->stats->expAA[EXP_CP9_GF][p]->dbsize, cm->stats->expAA[EXP_CP9_GF][p]->nrandhits, cm->stats->expAA[EXP_CP9_GF][p]->tailp);
	}
      /* currently either all exp tail stats are calc'ed or none */

      if (cm->flags & CMH_FILTER_STATS) /* FILTER stats are only possibly valid IF exp tail stats valid */
	{
	  fprintf(fp, "FT-LC    %d  %.5f  %d  %ld  %d\n", 
		  cm->stats->hfiA[FTHR_CM_LC]->ncut,  cm->stats->hfiA[FTHR_CM_LC]->F,
		  cm->stats->hfiA[FTHR_CM_LC]->N,     cm->stats->hfiA[FTHR_CM_LC]->dbsize,
		  cm->stats->hfiA[FTHR_CM_LC]->always_better_than_Smax);
	  fprintf(fp, "         ");
	  for(i = 0; i < cm->stats->hfiA[FTHR_CM_LC]->ncut; i++) fprintf(fp, "%10g ", cm->stats->hfiA[FTHR_CM_LC]->cm_E_cut[i]);
	  fprintf(fp, "\n");
	  fprintf(fp, "         ");
	  for(i = 0; i < cm->stats->hfiA[FTHR_CM_LC]->ncut; i++) fprintf(fp, "%10g ", cm->stats->hfiA[FTHR_CM_LC]->fwd_E_cut[i]);
	  fprintf(fp, "\n");

	  fprintf(fp, "FT-LI    %d  %.5f  %d  %ld  %d\n", 
		  cm->stats->hfiA[FTHR_CM_LI]->ncut,  cm->stats->hfiA[FTHR_CM_LI]->F,
		  cm->stats->hfiA[FTHR_CM_LI]->N,     cm->stats->hfiA[FTHR_CM_LI]->dbsize,
		  cm->stats->hfiA[FTHR_CM_LI]->always_better_than_Smax);
	  fprintf(fp, "         ");
	  for(i = 0; i < cm->stats->hfiA[FTHR_CM_LI]->ncut; i++) fprintf(fp, "%10g ", cm->stats->hfiA[FTHR_CM_LI]->cm_E_cut[i]);
	  fprintf(fp, "\n");
	  fprintf(fp, "         ");
	  for(i = 0; i < cm->stats->hfiA[FTHR_CM_LI]->ncut; i++) fprintf(fp, "%10g ", cm->stats->hfiA[FTHR_CM_LI]->fwd_E_cut[i]);
	  fprintf(fp, "\n");

	  fprintf(fp, "FT-GC    %d  %.5f  %d  %ld  %d\n", 
		  cm->stats->hfiA[FTHR_CM_GC]->ncut,  cm->stats->hfiA[FTHR_CM_GC]->F,
		  cm->stats->hfiA[FTHR_CM_GC]->N,     cm->stats->hfiA[FTHR_CM_GC]->dbsize,
		  cm->stats->hfiA[FTHR_CM_GC]->always_better_than_Smax);
	  fprintf(fp, "         ");
	  for(i = 0; i < cm->stats->hfiA[FTHR_CM_GC]->ncut; i++) fprintf(fp, "%10g ", cm->stats->hfiA[FTHR_CM_GC]->cm_E_cut[i]);
	  fprintf(fp, "\n");
	  fprintf(fp, "         ");
	  for(i = 0; i < cm->stats->hfiA[FTHR_CM_GC]->ncut; i++) fprintf(fp, "%10g ", cm->stats->hfiA[FTHR_CM_GC]->fwd_E_cut[i]);
	  fprintf(fp, "\n");

	  fprintf(fp, "FT-GI    %d  %.5f  %d  %ld  %d\n", 
		  cm->stats->hfiA[FTHR_CM_GI]->ncut,  cm->stats->hfiA[FTHR_CM_GI]->F,
		  cm->stats->hfiA[FTHR_CM_GI]->N,     cm->stats->hfiA[FTHR_CM_GI]->dbsize,
		  cm->stats->hfiA[FTHR_CM_GI]->always_better_than_Smax);
	  fprintf(fp, "         ");
	  for(i = 0; i < cm->stats->hfiA[FTHR_CM_GI]->ncut; i++) fprintf(fp, "%10g ", cm->stats->hfiA[FTHR_CM_GI]->cm_E_cut[i]);
	  fprintf(fp, "\n");
	  fprintf(fp, "         ");
	  for(i = 0; i < cm->stats->hfiA[FTHR_CM_GI]->ncut; i++) fprintf(fp, "%10g ", cm->stats->hfiA[FTHR_CM_GI]->fwd_E_cut[i]);
	  fprintf(fp, "\n");
	} /* currently either all filter threshold stats are calc'ed or none */
    }

  /* main model section */
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
	  for (x = 0; x < cm->abc->K; x++)
	    for (y = 0; y < cm->abc->K; y++)
	      fprintf(fp, "%6s ", prob2ascii(cm->e[v][x*cm->abc->K+y], cm->null[x]*cm->null[y]));
	}
      else if (cm->sttype[v] == ML_st || cm->sttype[v] == MR_st || cm->sttype[v] == IL_st || cm->sttype[v] == IR_st)
	{
	  for (x = 0; x < cm->abc->K; x++)
	    fprintf(fp, "%6s ", prob2ascii(cm->e[v][x], cm->null[x]));
	}
      fputs("\n", fp);
    }
  fputs("//\n", fp);
  return eslOK;
} 

static int  
read_ascii_cm(CMFILE *cmf, char *errbuf, ESL_ALPHABET **ret_abc, CM_t **ret_cm)
{
  int     status;
  CM_t   *cm;
  char   *buf;
  int     n;			/* length of buf */
  char   *s;
  int     M,N;			/* number of states, nodes in model */
  int     v,x,y,nd;		/* counters for states, events, nodes */
  char   *tok;
  int     toklen;
  int     exp_flags[EXP_NMODES]; /* keep track of which exp tails we've read */
  int     fthr_flags[FTHR_NMODES];/* keep track of which filter thresholds we've read */
  int     exp_mode;             /* index of exp tail info               */
  int     fthr_mode;            /* HMM filter threshold info       */
  int     have_exps;            /* for checking we get 0 or all exp tails*/
  int     have_fthrs;           /* for checking we get 0 or all fthrs */
  int     have_ga = FALSE;      /* we have GA cutoff, needed b/c we can't set cm->flags until after CreateCMBody() call */
  int     have_tc = FALSE;      /* we have TC cutoff, needed b/c we can't set cm->flags until after CreateCMBody() call */
  int     have_nc = FALSE;      /* we have NC cutoff, needed b/c we can't set cm->flags until after CreateCMBody() call */
  int     p;                    /* counter for partitions          */
  int     gc;                   /* counter over gc contents        */
  int     i;                    /* counter over exp_modes for exp tails */
  int     alphabet_type;        /* type of ESL_ALPHABET */
  ESL_ALPHABET *abc = NULL;
  int     read_nstates = FALSE; /* TRUE once we've read the number of states */
  int     read_nnodes  = FALSE; /* TRUE once we've read the number of nodes */
  int     read_clen = FALSE;
  int     clen = 0;

  cm  = NULL;
  buf = NULL;
  n   = 0;
  for(i = 0; i < EXP_NMODES; i++)  exp_flags[i] = FALSE;
  for(i = 0; i < FTHR_NMODES; i++) fthr_flags[i] = FALSE;

  if (feof(cmf->f) || esl_fgets(&buf, &n, cmf->f) != eslOK) { /* end of file, free buf and return eslEOF */
    if(buf != NULL) free(buf);
    return eslEOF;
  }
  
  if (strncmp(buf, "INFERNAL-1", 10) != 0)                 goto FAILURE;

  /* Parse the header information
   * These are all tag/value. 
   * Ignore unknown tags (forward compatibility). 
   */
  cm = CreateCMShell();
  M  = N = -1;
  while (esl_fgets(&buf, &n, cmf->f) != eslEOF) 
    {
      s   = buf;
      if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;
      else if (strcmp(tok, "NAME") == 0) 
	{
	  if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;
	  if ((esl_strdup(tok, toklen, &(cm->name)))   != eslOK) goto ERROR;
	}
      else if (strcmp(tok, "ACC") == 0) 
	{
	  if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;
	  if ((esl_strdup(tok, toklen, &(cm->acc)))    != eslOK) goto ERROR;
	}
      else if (strcmp(tok, "DESC") == 0) 
	{
	  if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;
	  esl_strdup(tok, toklen, &(cm->desc));
	}
      else if (strcmp(tok, "GA") == 0) 
	{
	  if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;
	  cm->ga = atof(tok);
	  have_ga = TRUE;
	}
      else if (strcmp(tok, "TC") == 0) 
	{
	  if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;
	  cm->tc = atof(tok);
	  have_tc = TRUE;
	}
      else if (strcmp(tok, "NC") == 0) 
	{
	  if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;
	  cm->nc = atof(tok);
	  have_nc = TRUE;
	}
      else if (strcmp(tok, "STATES") == 0) 
	{
	  if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;
	  M = atoi(tok);
	  read_nstates = TRUE;
	}
      else if (strcmp(tok, "NODES") == 0) 
	{
	  if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;
	  N = atoi(tok);
	  read_nnodes = TRUE;
	}
      else if (strcmp(tok, "ALPHABET") == 0) 
	{
	  if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;
	  alphabet_type = atoi(tok);
	  /* Set or verify alphabet. */
	  if (*ret_abc == NULL)	{	/* still unknown: set it, pass control of it back to caller */
	    if ((abc = esl_alphabet_Create(alphabet_type)) == NULL)       { status = eslEMEM;      goto FAILURE; }
	  } else {			/* already known: check it */
	    abc = *ret_abc;
	    if ((*ret_abc)->type != alphabet_type)                        { status = eslEINCOMPAT; goto FAILURE; }
	  }
	  /* Now we have the alphabet and we should have N and M, so we can build the
	   * full model, and set the alphabet (which we need to do before alloc'ing/setting
	   * the null model */
	  if(! (read_nstates && read_nnodes))
	    {
	      printf("ERROR, STATES and NODES lines should precede alphabet line");
	      goto FAILURE;
	    }
	  CreateCMBody(cm, N, M, abc);
	}	    
      else if (strcmp(tok, "ELSELF") == 0) 
	{
	  if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;
	  cm->el_selfsc = atof(tok);
	}
      else if (strcmp(tok, "WBETA") == 0) 
	{
	  if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;
	  cm->beta_W = (double) atof(tok);
	  cm->beta_qdb = cm->beta_W;
	}
      else if (strcmp(tok, "NSEQ") == 0) 
	{
	  if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;
	  cm->nseq = atoi(tok);
	}
      else if (strcmp(tok, "EFFNSEQ") == 0) 
	{
	  if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;
	  cm->eff_nseq = atof(tok);
	}
      else if (strcmp(tok, "CLEN") == 0) 
	{
	  if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;
	  clen = atoi(tok); /* we'll compare this to what we calculate at end of func */
	  read_clen = TRUE;
	}
      /* comlog info, careful, we want the full line, so a token becomes a full line */
      else if (strcmp(tok, "BCOM") == 0) 
	{
	  while(isspace((int) (*s))) s++; /* chew up leading whitespace */
	  if ((esl_strtok(&s, "\n", &tok, &toklen)) != eslOK) goto FAILURE;
	  if(cm->comlog->bcom != NULL) free(cm->comlog->bcom);
	  esl_strdup(tok, toklen, &(cm->comlog->bcom));
	}
      else if (strcmp(tok, "BDATE") == 0) 
	{
	  while(isspace((int) (*s))) s++; /* chew up leading whitespace */
	  if ((esl_strtok(&s, "\n", &tok, &toklen)) != eslOK) goto FAILURE;
	  if(cm->comlog->bdate != NULL) free(cm->comlog->bdate);
	  esl_strdup(tok, toklen, &(cm->comlog->bdate));
	}
      else if (strcmp(tok, "CCOM") == 0) 
	{
	  while(isspace((int) (*s))) s++; /* chew up leading whitespace */
	  if ((esl_strtok(&s, "\n", &tok, &toklen)) != eslOK) goto FAILURE;
	  if(cm->comlog->ccom != NULL) free(cm->comlog->ccom);
	  esl_strdup(tok, toklen, &(cm->comlog->ccom));
	}
      else if (strcmp(tok, "CDATE") == 0) 
	{
	  while(isspace((int) (*s))) s++; /* chew up leading whitespace */
	  if ((esl_strtok(&s, "\n", &tok, &toklen)) != eslOK) goto FAILURE;
	  if(cm->comlog->cdate != NULL) free(cm->comlog->cdate);
	  esl_strdup(tok, toklen, &(cm->comlog->cdate));
	}
      else if (strcmp(tok, "NULL") == 0) 
	{
	  if(cm->abc == NULL) goto FAILURE;
	  /* cm-> null already allocated in CreateCMBody() */
	  for (x = 0; x < abc->K; x++)
	    {
	      if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;
	      cm->null[x] = ascii2prob(tok, (1./(float) abc->K));
	    }
	}
      /* exp tail distribution information */
      else if (strcmp(tok, "PART") == 0) 
	{
	  if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;
	  if (! is_integer(tok))                                      goto FAILURE;
	  /* First token is num partitions, allocate cmstats object based on this */
	  cm->stats = AllocCMStats(atoi(tok));
	  for(p = 0; p < cm->stats->np; p++)
	    {
	      /* there are 2 * cm->stats->np tokens left on this line,
	       * (ps, pe) pairs for each partition */
	      if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;
	      if (! is_integer(tok))                                      goto FAILURE;
	      cm->stats->ps[p] = atoi(tok);
	      if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;
	      if (! is_integer(tok))                                      goto FAILURE;
	      cm->stats->pe[p] = atoi(tok);
	    }
	  /* Now set the gc2p GC content to partition map, 
	   * [0..GC_SEGMENTS], telling which partition each belongs to */
	  gc = 0;
	  for(p = 0; p < cm->stats->np; p++)
	    {
	      if(cm->stats->ps[p] != gc)                   goto FAILURE;
	      while(gc <= cm->stats->pe[p]) 
		cm->stats->gc2p[gc++] = p;
	    }
	  if(gc != GC_SEGMENTS)                         goto FAILURE;
	}
      /* exp tail info */
      else if (strncmp(tok, "E-", 2) == 0) 
      {				
	/* determine which exp tail we're reading */
	if      (strncmp(tok+2, "LC", 2) == 0) 
	  exp_mode = EXP_CM_LC;
	else if (strncmp(tok+2, "GC", 2) == 0) 
	  exp_mode = EXP_CM_GC;
	else if (strncmp(tok+2, "LI", 2) == 0) 
	  exp_mode = EXP_CM_LI;
	else if (strncmp(tok+2, "GI", 2) == 0) 
	  exp_mode = EXP_CM_GI;
	else if (strncmp(tok+2, "LV", 2) == 0) 
	  exp_mode = EXP_CP9_LV;
	else if (strncmp(tok+2, "GV", 2) == 0) 
	  exp_mode = EXP_CP9_GV;
	else if (strncmp(tok+2, "LF", 2) == 0) 
	  exp_mode = EXP_CP9_LF;
	else if (strncmp(tok+2, "GF", 2) == 0) 
	  exp_mode = EXP_CP9_GF;
	else                                         goto FAILURE;

	/* now we know what exp tail we're reading, read it */
	if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;
	if (! is_integer(tok))                                 goto FAILURE;
	p = atoi(tok);
	if (p >= cm->stats->np)                                goto FAILURE;

	if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;
	if (! is_real(tok))                                    goto FAILURE;
	cm->stats->expAA[exp_mode][p]->lambda = atof(tok);

	if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;
	if (! is_real(tok))                                    goto FAILURE;
	cm->stats->expAA[exp_mode][p]->mu_extrap = atof(tok);

	if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;
	if (! is_real(tok))                                    goto FAILURE;
	cm->stats->expAA[exp_mode][p]->mu_orig = atof(tok);

	if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;
	if (! is_integer(tok))                                 goto FAILURE;
	cm->stats->expAA[exp_mode][p]->dbsize = (long) atoi(tok);

	if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;
	if (! is_integer(tok))                                 goto FAILURE;
	cm->stats->expAA[exp_mode][p]->nrandhits = atoi(tok);

	if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;
	if (! is_real(tok))                                    goto FAILURE;
	cm->stats->expAA[exp_mode][p]->tailp = atof(tok);

	cm->stats->expAA[exp_mode][p]->cur_eff_dbsize = (long) (cm->stats->expAA[exp_mode][p]->nrandhits);
	/* Previous line is to set cur_eff_dbsize as if database was of size cm->stats->expAA[p]->dbsize, we 
	 * act as if the max hits we'll see is nrandhits, the number of hits we saw in cmcalibrate,
	 * so this is the highest possible E-value we can get.
	 * cur_eff_dbsize will be updated in cmsearch for whatever the target database size is. */
	cm->stats->expAA[exp_mode][p]->is_valid = TRUE; /* set valid flag */
	exp_flags[exp_mode] = TRUE;
      }
      /* best filter threshold info */
      else if (strncmp(tok, "FT-", 3) == 0) 
      {				
	/* cm->stats should've been alloc'ed when exp tails were read */
	if(cm->stats == NULL) 	                     goto FAILURE;
	/* determine which filter threshold we're reading */
	if (strncmp(tok+3, "LC", 2) == 0) 
	  fthr_mode = FTHR_CM_LC;
	else if (strncmp(tok+3, "GC", 2) == 0) 
	  fthr_mode = FTHR_CM_GC;
	else if (strncmp(tok+3, "LI", 2) == 0) 
	  fthr_mode = FTHR_CM_LI;
	else if (strncmp(tok+3, "GI", 2) == 0) 
	  fthr_mode = FTHR_CM_GI;
	else                                         goto FAILURE;

	/* now we know what mode we're reading, read it */
	if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;
	if (! is_integer(tok))                                 goto FAILURE;
	cm->stats->hfiA[fthr_mode]->ncut = atoi(tok);
	if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;
	if (! is_real(tok))                                    goto FAILURE;
	cm->stats->hfiA[fthr_mode]->F = atof(tok);
	if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;
	if (! is_integer(tok))                                 goto FAILURE;
	cm->stats->hfiA[fthr_mode]->N = atoi(tok);
	if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;
	if (! is_integer(tok))                                 goto FAILURE;
	cm->stats->hfiA[fthr_mode]->dbsize = (long) atoi(tok);
	if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;
	if (! is_real(tok))                                    goto FAILURE;
	cm->stats->hfiA[fthr_mode]->always_better_than_Smax = atoi(tok);

	/* alloc for, and read a new line, the CM cut points */
	ESL_ALLOC(cm->stats->hfiA[fthr_mode]->cm_E_cut,  sizeof(float) * cm->stats->hfiA[fthr_mode]->ncut);
	if (esl_fgets(&buf, &n, cmf->f) != eslOK) goto FAILURE;
	s = buf;
	for(i = 0; i < cm->stats->hfiA[fthr_mode]->ncut; i++) { 
	  if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;
	  if (! is_real(tok))                                    goto FAILURE;
	  cm->stats->hfiA[fthr_mode]->cm_E_cut[i] = atof(tok);
	}

	ESL_ALLOC(cm->stats->hfiA[fthr_mode]->fwd_E_cut, sizeof(float) * cm->stats->hfiA[fthr_mode]->ncut);
	if (esl_fgets(&buf, &n, cmf->f) != eslOK) goto FAILURE;
	s = buf;
	for(i = 0; i < cm->stats->hfiA[fthr_mode]->ncut; i++) { 
	  if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;
	  if (! is_real(tok))                                    goto FAILURE;
	  cm->stats->hfiA[fthr_mode]->fwd_E_cut[i] = atof(tok);
	}
	cm->stats->hfiA[fthr_mode]->is_valid = TRUE;
	fthr_flags[fthr_mode] = TRUE;
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

  /* if we have any exp tail stats, we (currently) require all of them */
  have_exps = exp_flags[0];
  for(exp_mode = 1; exp_mode < EXP_NMODES; exp_mode++)
    if(((have_exps && (!exp_flags[exp_mode]))) ||
       ((!have_exps) && (exp_flags[exp_mode])))
      goto FAILURE;
  
  /* if we have any filter stats, we (currently) require all of them */
  have_fthrs = fthr_flags[0];
  for(i = 0; i < FTHR_NMODES; i++) {
    if(((have_fthrs && (!fthr_flags[i]))) || ((!have_fthrs) && (fthr_flags[i]))) goto FAILURE;
  }
  /* if we have exp tail stats we must have filter thresholds stats, 
   * and if we have filter threshold stats we must have exp tail stats.
   */
  if(have_exps  && !have_fthrs) goto FAILURE;
  if(!have_exps &&  have_fthrs) goto FAILURE;

  /* Main model section. 
   */
  CMZero(cm);
  if(have_exps)  cm->flags |= CMH_EXPTAIL_STATS;
  if(have_fthrs) cm->flags |= CMH_FILTER_STATS;
  if(have_ga)    cm->flags |= CMH_GA;
  if(have_tc)    cm->flags |= CMH_TC;
  if(have_nc)    cm->flags |= CMH_NC;
  nd = -1;
  cm->clen = 0;
  for (v = 0; v < cm->M; v++)
    {
      if (esl_fgets(&buf, &n, cmf->f) != eslOK) goto FAILURE;
      s = buf;
      if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;      
      
      /* Ah, a node line. Process it and get the following line.
       */
      if (*tok == '[') 
	{
	  if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;      
	  if ((x = NodeCode(tok)) == -1)                        goto FAILURE;
	  if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;      
	  if (!is_integer(tok))                                      goto FAILURE;
	  nd = atoi(tok);
	  cm->ndtype[nd]  = x;
	  if(cm->ndtype[nd] == MATP_nd) cm->clen+=2;
	  else if(cm->ndtype[nd] == MATL_nd) cm->clen++;
	  else if(cm->ndtype[nd] == MATR_nd) cm->clen++;
	  cm->nodemap[nd] = v;

	  if (esl_fgets(&buf, &n, cmf->f) != eslOK)              goto FAILURE;
	  s = buf;
	  if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;      
	}

      /* Process state line.
       */
      cm->sttype[v] = StateCode(tok);
      if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;      
      if (! is_integer(tok))                                 goto FAILURE;
      if (atoi(tok) != v)                                    goto FAILURE;
      if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;      
      if (! is_integer(tok))                                 goto FAILURE;
      cm->plast[v] = atoi(tok);
      if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;      
      if (! is_integer(tok))                                 goto FAILURE;
      cm->pnum[v] = atoi(tok);
      if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;      
      if (! is_integer(tok))                                 goto FAILURE;
      cm->cfirst[v] = atoi(tok);
      if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;      
      if (! is_integer(tok))                                 goto FAILURE;
      cm->cnum[v] = atoi(tok);
				/* Transition probabilities. */
      if (cm->sttype[v] != B_st) 
	{
	  for (x = 0; x < cm->cnum[v]; x++)
	    {
	      if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;      
	      if (! is_real(tok) && *tok != '*')                      goto FAILURE;
	      cm->t[v][x] = ascii2prob(tok, 1.);
	    }
	}
				/* Emission probabilities. */
      if (cm->sttype[v] == ML_st || cm->sttype[v] == MR_st ||
	  cm->sttype[v] == IL_st || cm->sttype[v] == IR_st)
	{
	  for (x = 0; x < cm->abc->K; x++)
	    {
	      if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;      
	      if (! is_real(tok) && *tok != '*')                     goto FAILURE;
	      cm->e[v][x] = ascii2prob(tok, cm->null[x]);
	    }
	}
      else if (cm->sttype[v] == MP_st) 
	{
	  for (x = 0; x < cm->abc->K; x++)
	    for (y = 0; y < cm->abc->K; y++)
	      {
		if ((esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto FAILURE;      
		if (! is_real(tok) && *tok != '*')                     goto FAILURE;
		cm->e[v][x*cm->abc->K+y] = ascii2prob(tok, cm->null[x]*cm->null[y]);
	      }
	} 

      cm->ndidx[v] = nd;
      cm->stid[v]  = DeriveUniqueStateCode(cm->ndtype[nd], cm->sttype[v]);
    } /* end of loop over states */

  /* Advance to record separator
   */
  while (esl_fgets(&buf, &n, cmf->f) != eslEOF) 
    if (strncmp(buf, "//", 2) == 0) 
      break;

  /* EPN 10.29.06 Remove the sole source of CM ambiguities. Find and detach insert states
   *              that are 1 state before an END_E.  */
  cm_find_and_detach_dual_inserts(cm, 
				  FALSE, /* Don't check END_E-1 states have 0 counts, they may not if 
					  * an old version (0.7 or earlier) of cmbuild was used, or  
					  * cmbuild --nodetach  was used to build the CM  */
				  TRUE); /* Detach the states by setting trans probs into them as 0.0   */

  /* check that the clen we calc'ed is the same as the CLEN line said */
  if (read_clen && clen != cm->clen) 
    {
      printf("ERROR, calculated consensus length %d does not equal read CLEN: %d.\n", cm->clen, clen);
      goto FAILURE;
    }

  /* Success.
   * Renormalize the CM, and return.
   */
  CMRenormalize(cm);

  if (buf != NULL) free(buf);
  if (*ret_abc == NULL) *ret_abc = abc;	/* pass our new alphabet back to caller, if caller didn't know it already */
  *ret_cm = cm;
  return eslOK;

 FAILURE:
  if (buf != NULL) free(buf);
  *ret_cm = NULL;
  ESL_FAIL(eslEFORMAT, errbuf, "Error reading the cmfile. Is it corrupt or built with a pre-1.0 cmbuild?"); 
  return status; /* NEVERREACHED */

 ERROR:
  ESL_FAIL(eslEMEM, errbuf, "Error ran out of memory reading the cmfile.");
  return status; /* NEVERREACHED */
}

/* Function: write_binary_cm()
 * Date:     SRE, Thu Aug  3 12:05:30 2000 [St. Louis]
 *
 * Purpose:  Write a CM in binary format.
 *
 */
static int
write_binary_cm(FILE *fp, CM_t *cm, char *errbuf)
{
  int v, i ,p;
  int has_exp, has_fthr;
  int has_ga, has_tc, has_nc;
  int atype;
  atype = cm->abc->type;

  fwrite((char *) &(v01magic), sizeof(unsigned int), 1, fp);

  /* These have to go first, so we know how big of a CM
   * to allocate when we go to read a file.
   */
  tagged_fwrite(CMIO_M,            &cm->M,          sizeof(int),   1, fp);
  tagged_fwrite(CMIO_NODES,        &cm->nodes,      sizeof(int),   1, fp);  
  tagged_fwrite(CMIO_ALPHABETTYPE, &atype,          sizeof(int),   1, fp);

  tagged_bin_string_write(CMIO_NAME, cm->name,  fp);
  tagged_bin_string_write(CMIO_ACC,  cm->acc,   fp);
  tagged_bin_string_write(CMIO_DESC, cm->desc,  fp);
    /* Rfam cutoffs */
  if (!(cm->flags & CMH_GA))  { 
    has_ga = FALSE;
    tagged_fwrite(CMIO_HASGA, &has_ga, sizeof(int),  1, fp);  /* put a 0 to indicate no GA cutoff */
  }
  else {
    has_ga = TRUE;
    tagged_fwrite(CMIO_HASGA, &has_ga, sizeof(int),  1, fp);  /* put a 1 to indicate GA cutoff is next */
    tagged_fwrite(CMIO_GA,    &cm->ga, sizeof(int),  1, fp);
  }    
  if (!(cm->flags & CMH_TC))  { 
    has_tc = FALSE;
    tagged_fwrite(CMIO_HASTC, &has_tc, sizeof(int),  1, fp);  /* put a 0 to indicate no TC cutoff */
  }
  else {
    has_tc = TRUE;
    tagged_fwrite(CMIO_HASTC, &has_tc, sizeof(int),  1, fp);  /* put a 1 to indicate TC cutoff is next */
    tagged_fwrite(CMIO_TC,    &cm->tc, sizeof(int),  1, fp);
  }    
  if (!(cm->flags & CMH_NC))  { 
    has_nc = FALSE;
    tagged_fwrite(CMIO_HASNC, &has_nc, sizeof(int),  1, fp);  /* put a 0 to indicate no NC cutoff */
  }
  else {
    has_nc = TRUE;
    tagged_fwrite(CMIO_HASNC, &has_nc, sizeof(int),  1, fp);  /* put a 1 to indicate NC cutoff is next */
    tagged_fwrite(CMIO_NC,    &cm->nc, sizeof(int),  1, fp);
  }    

  tagged_fwrite(CMIO_ELSELFSC,    &cm->el_selfsc,   sizeof(float), 1, fp);  
  tagged_fwrite(CMIO_WBETA,       &cm->beta_W,      sizeof(double),1, fp);  
  tagged_fwrite(CMIO_NSEQ,        &cm->nseq,        sizeof(int),   1, fp);  
  tagged_fwrite(CMIO_EFFNSEQ,     &cm->eff_nseq,    sizeof(float), 1, fp);  
  tagged_fwrite(CMIO_CLEN,        &cm->clen,        sizeof(int),   1, fp);  

  /* cm->comlog, the creation dates and command lines used to build/calibrate the model */
  tagged_bin_string_write(CMIO_BCOM,   cm->comlog->bcom,  fp);
  tagged_bin_string_write(CMIO_BDATE,  cm->comlog->bdate, fp);
  tagged_bin_string_write(CMIO_CCOM,   cm->comlog->ccom, fp);
  tagged_bin_string_write(CMIO_CDATE,  cm->comlog->cdate,fp);
  /* null, background distro */
  tagged_fwrite(CMIO_NULL,         cm->null,       sizeof(float), cm->abc->K, fp);

  /* exp tail stats */
  if (!(cm->flags & CMH_EXPTAIL_STATS))
    {
      has_exp = FALSE;
      tagged_fwrite(CMIO_HASEXP,     &has_exp,  sizeof(int),  1, fp);  /* put a 0 to indicate no exp tail stats */
    }
  else /* (cm->flags & CMH_EXPTAILL_STATS), if this flag is up, ALL exp tail stats are valid */
    {
      has_exp = TRUE;
      tagged_fwrite(CMIO_HASEXP,  &has_exp,         sizeof(int),  1, fp);  /* put a 1 to indicate valid exp tail stats */
      tagged_fwrite(CMIO_NPART,   &cm->stats->np,   sizeof(int),  1, fp);  
      tagged_fwrite(CMIO_PARTS,   cm->stats->ps,    sizeof(int),  cm->stats->np, fp);  
      tagged_fwrite(CMIO_PARTE,   cm->stats->pe,    sizeof(int),  cm->stats->np, fp);  
      for(i = 0; i < EXP_NMODES; i++)
	{
	  for(p = 0; p < cm->stats->np; p++)
	    {
	      tagged_fwrite(CMIO_EXPLAMBDA, &cm->stats->expAA[i][p]->lambda,    sizeof(double),  1, fp);
	      tagged_fwrite(CMIO_EXPMUE,    &cm->stats->expAA[i][p]->mu_extrap, sizeof(double),  1, fp);
	      tagged_fwrite(CMIO_EXPMUO,    &cm->stats->expAA[i][p]->mu_orig,   sizeof(double),  1, fp);
	      tagged_fwrite(CMIO_EXPDBSIZE, &cm->stats->expAA[i][p]->dbsize,    sizeof(long),    1, fp);
	      tagged_fwrite(CMIO_EXPNHITS,  &cm->stats->expAA[i][p]->nrandhits, sizeof(int),     1, fp);
	      tagged_fwrite(CMIO_EXPTAILP,  &cm->stats->expAA[i][p]->tailp,     sizeof(double),  1, fp);
	    }
	}
    }
  /* HMM filter threshold stats */
  if (!(cm->flags & CMH_FILTER_STATS))
    {
      has_fthr = FALSE;
      tagged_fwrite(CMIO_HASFILTER,  &has_fthr,   sizeof(int),  1, fp);  /* put a 0 to indicate no HMM filter stats */
    } 
  else /* (cm->flags & CMH_FILTERSTATS), check to make sure exp tail stats are also valid, they should be */
    {
      has_fthr = TRUE;
      if(! (cm->flags & CMH_EXPTAIL_STATS)) cm_Fail("writing binary CM file, filter stats were valid, but exp tail stats were not, this shouldn't happen.");
      tagged_fwrite(CMIO_HASFILTER,  &has_fthr,   sizeof(int),  1, fp);  /* put a 1 to indicate valid HMM filter stats */
      for(i = 0; i < FTHR_NMODES; i++)
	{
	  tagged_fwrite(CMIO_FTHRNCUT,   &cm->stats->hfiA[i]->ncut,                    sizeof(int),   1, fp);      
	  tagged_fwrite(CMIO_FTHRF,      &cm->stats->hfiA[i]->F,                       sizeof(float), 1, fp);      
	  tagged_fwrite(CMIO_FTHRN,      &cm->stats->hfiA[i]->N,                       sizeof(int),   1, fp);      
	  tagged_fwrite(CMIO_FTHRDB,     &cm->stats->hfiA[i]->dbsize,                  sizeof(long),  1, fp);      
	  tagged_fwrite(CMIO_FTHRABTS,   &cm->stats->hfiA[i]->always_better_than_Smax, sizeof(int),   1, fp);      
	  tagged_fwrite(CMIO_FTHRCMECUT,  cm->stats->hfiA[i]->cm_E_cut,                sizeof(float), cm->stats->hfiA[i]->ncut, fp);      
	  tagged_fwrite(CMIO_FTHRFWDECUT, cm->stats->hfiA[i]->fwd_E_cut,               sizeof(float), cm->stats->hfiA[i]->ncut, fp);      
	}
    }

  /* main model section */
  tagged_fwrite(CMIO_STTYPE,       cm->sttype,     sizeof(char),  cm->M, fp);
  tagged_fwrite(CMIO_NDIDX,        cm->ndidx,      sizeof(int),   cm->M, fp);  
  tagged_fwrite(CMIO_STID,         cm->stid,       sizeof(char),  cm->M, fp);  
  tagged_fwrite(CMIO_CFIRST,       cm->cfirst,     sizeof(int),   cm->M, fp); 
  tagged_fwrite(CMIO_CNUM,         cm->cnum,       sizeof(int),   cm->M, fp);  
  tagged_fwrite(CMIO_PLAST,        cm->plast,      sizeof(int),   cm->M, fp); 
  tagged_fwrite(CMIO_PNUM,         cm->pnum,       sizeof(int),   cm->M, fp);  
  tagged_fwrite(CMIO_NODEMAP,      cm->nodemap,    sizeof(int),   cm->nodes, fp);
  tagged_fwrite(CMIO_NDTYPE,       cm->ndtype,     sizeof(char),  cm->nodes, fp);
  for (v = 0; v < cm->M; v++) {
    tagged_fwrite(CMIO_T, cm->t[v], sizeof(float), MAXCONNECT, fp);
    tagged_fwrite(CMIO_E, cm->e[v], sizeof(float), cm->abc->K*cm->abc->K, fp);
  }
  tagged_fwrite(CMIO_END_DATA, NULL, 0, 0, fp);

  /* Note: begin, end, and flags not written out. Local alignment is
   * run-time configuration right now.
   */
  return eslOK;
}


/* Function: read_binary_cm()
 * Date:     SRE, Thu Aug  3 13:39:09 2000 [St. Louis]
 *
 * Purpose:  Read a CM from disk.
 */
static int
read_binary_cm(CMFILE *cmf, char *errbuf, ESL_ALPHABET **ret_abc, CM_t **ret_cm)
{
  FILE         *fp;
  CM_t         *cm;
  unsigned int  magic;
  int           M;
  int           nodes;
  int           alphabet_type;
  int           v;
  int           has_exp;
  int           has_fthr;
  int           has_ga, has_tc, has_nc;
  int           np;
  int           i, p, gc;
  ESL_ALPHABET *abc = NULL;
  int           status;

  cm = NULL;
  fp = cmf->f;
  if (feof(fp)) return eslEOF;
  if (! fread((char *) &magic, sizeof(unsigned int), 1, fp)) return eslEOF;
  if (magic != v01magic) goto FAILURE;
  
  if (! tagged_fread(CMIO_M,     (void *) &M,     sizeof(int), 1, fp)) goto FAILURE;
  if (! tagged_fread(CMIO_NODES, (void *) &nodes, sizeof(int), 1, fp)) goto FAILURE;
  if (! tagged_fread(CMIO_ALPHABETTYPE,(void *) &alphabet_type, sizeof(int), 1, fp)) goto FAILURE;

  /* Set or verify alphabet. */
  if (*ret_abc == NULL)	{	/* still unknown: set it, pass control of it back to caller */
    if ((abc = esl_alphabet_Create(alphabet_type)) == NULL)       { status = eslEMEM;      goto FAILURE; }
  } else {			/* already known: check it */
    abc = *ret_abc;
    if ((*ret_abc)->type != alphabet_type)                        { status = eslEINCOMPAT; goto FAILURE; }
  }
  cm = CreateCM(nodes, M, abc);

  if (! tagged_bin_string_read(CMIO_NAME, &(cm->name),  fp)) goto FAILURE;
  if (! tagged_bin_string_read(CMIO_ACC,  &(cm->acc),   fp)) goto FAILURE;
  if (! tagged_bin_string_read(CMIO_DESC, &(cm->desc),  fp)) goto FAILURE;
  /* we might have any combo of Rfam cutoffs */
  if (! tagged_fread(CMIO_HASGA,    (void *) &(has_ga),     sizeof(int),   1,         fp))    goto FAILURE;
  if(has_ga) { 
    if (! tagged_fread(CMIO_GA,    (void *) &(cm->ga),       sizeof(int),   1,         fp))    goto FAILURE;
    cm->flags |= CMH_GA;
  }
  if (! tagged_fread(CMIO_HASTC,    (void *) &(has_tc),     sizeof(int),   1,         fp))    goto FAILURE;
  if(has_tc) { 
    if (! tagged_fread(CMIO_TC,    (void *) &(cm->tc),       sizeof(int),   1,         fp))    goto FAILURE;
    cm->flags |= CMH_TC;
  }
  if (! tagged_fread(CMIO_HASNC,    (void *) &(has_nc),     sizeof(int),   1,         fp))    goto FAILURE;
  if(has_nc) { 
    if (! tagged_fread(CMIO_NC,    (void *) &(cm->nc),       sizeof(int),   1,         fp))    goto FAILURE;
    cm->flags |= CMH_NC;
  }
  if (! tagged_fread(CMIO_ELSELFSC, (void *) &(cm->el_selfsc), sizeof(float), 1, fp)) goto FAILURE;
  if (! tagged_fread(CMIO_WBETA,    (void *) &(cm->beta_W),    sizeof(double),1, fp)) goto FAILURE;
  cm->beta_qdb = cm->beta_W;
  if (! tagged_fread(CMIO_NSEQ,     (void *) &(cm->nseq),      sizeof(int),   1, fp)) goto FAILURE;
  if (! tagged_fread(CMIO_EFFNSEQ,  (void *) &(cm->eff_nseq),  sizeof(float), 1, fp)) goto FAILURE;
  if (! tagged_fread(CMIO_CLEN,     (void *) &(cm->clen),      sizeof(int),   1, fp)) goto FAILURE;

  /* comlog info */
  if (! tagged_bin_string_read(CMIO_BCOM,   &(cm->comlog->bcom),   fp)) goto FAILURE;
  if (! tagged_bin_string_read(CMIO_BDATE,  &(cm->comlog->bdate),  fp)) goto FAILURE;
  if (! tagged_bin_string_read(CMIO_CCOM,   &(cm->comlog->ccom),  fp)) goto FAILURE;
  if (! tagged_bin_string_read(CMIO_CDATE,  &(cm->comlog->cdate), fp)) goto FAILURE;

  /* null distro */
  if (! tagged_fread(CMIO_NULL,         (void *) cm->null,       sizeof(float), cm->abc->K, fp))    goto FAILURE;
  /* We might have exp tail stats */
  if (! tagged_fread(CMIO_HASEXP,       (void *) &(has_exp),     sizeof(int),   1,         fp))    goto FAILURE;
  if(has_exp)
    {
      /* First is num partitions, allocate cmstats object based on this */
      if (! tagged_fread(CMIO_NPART,     (void *) &(np),         sizeof(int),       1,        fp))     goto FAILURE;
      cm->stats = AllocCMStats(np);
      if (! tagged_fread(CMIO_PARTS,     (void *) cm->stats->ps, sizeof(int),      np,        fp))     goto FAILURE;
      if (! tagged_fread(CMIO_PARTE,     (void *) cm->stats->pe, sizeof(int),      np,        fp))     goto FAILURE;
      /* Now set the gc2p GC content to partition map, 
       * [0..GC_SEGMENTS], telling which partition each belongs to */
      gc = 0;
      for(p = 0; p < cm->stats->np; p++) {
	if(cm->stats->ps[p] != gc)                   goto FAILURE;
	while(gc <= cm->stats->pe[p]) 
	  cm->stats->gc2p[gc++] = p;
      }
      if(gc != GC_SEGMENTS)                         goto FAILURE;

      for(i = 0; i < EXP_NMODES; i++)
	{
	  for(p = 0; p < cm->stats->np; p++)
	    {
	      if (! tagged_fread(CMIO_EXPLAMBDA, (void *) &(cm->stats->expAA[i][p]->lambda),   sizeof(double), 1, fp)) goto FAILURE;
	      if (! tagged_fread(CMIO_EXPMUE,    (void *) &(cm->stats->expAA[i][p]->mu_extrap),sizeof(double), 1, fp)) goto FAILURE;
	      if (! tagged_fread(CMIO_EXPMUO,    (void *) &(cm->stats->expAA[i][p]->mu_orig),  sizeof(double), 1, fp)) goto FAILURE;
	      if (! tagged_fread(CMIO_EXPDBSIZE, (void *) &(cm->stats->expAA[i][p]->dbsize),   sizeof(long),   1, fp)) goto FAILURE;
	      if (! tagged_fread(CMIO_EXPNHITS,  (void *) &(cm->stats->expAA[i][p]->nrandhits),sizeof(int),    1, fp)) goto FAILURE;
	      if (! tagged_fread(CMIO_EXPTAILP,  (void *) &(cm->stats->expAA[i][p]->tailp),    sizeof(double), 1, fp)) goto FAILURE;
	      cm->stats->expAA[i][p]->cur_eff_dbsize = (long) (cm->stats->expAA[i][p]->nrandhits);
	      /* Previous line is to set cur_eff_dbsize as if database was of size cm->stats->expAA[p]->dbsize, we 
	       * act as if the max hits we'll see is nrandhits, the number of hits we saw in cmcalibrate,
	       * so this is the highest possible E-value we can get.
	       * cur_eff_dbsize will be updated in cmsearch for whatever the target database size is. */
	      cm->stats->expAA[i][p]->is_valid = TRUE; /* set valid flag */
	    }
	}
      cm->flags |= CMH_EXPTAIL_STATS;
    }
  if (! tagged_fread(CMIO_HASFILTER,     (void *) &(has_fthr), sizeof(int),         1,        fp))     goto FAILURE;
  /* We might have HMM filter threshold stats */
  if(has_fthr)
    {
      if(! has_exp) goto FAILURE; /* filter threshold stats should only exist if exp tail stats exist */
      for(i = 0; i < FTHR_NMODES; i++)
	{
	  if (! tagged_fread(CMIO_FTHRNCUT,    (void *) &(cm->stats->hfiA[i]->ncut),                    sizeof(int),   1, fp)) goto FAILURE;
	  if (! tagged_fread(CMIO_FTHRF,       (void *) &(cm->stats->hfiA[i]->F),                       sizeof(float), 1, fp)) goto FAILURE;
	  if (! tagged_fread(CMIO_FTHRN,       (void *) &(cm->stats->hfiA[i]->N),                       sizeof(int),   1, fp)) goto FAILURE;
	  if (! tagged_fread(CMIO_FTHRDB,      (void *) &(cm->stats->hfiA[i]->dbsize),                  sizeof(long),  1, fp)) goto FAILURE;
	  if (! tagged_fread(CMIO_FTHRABTS,    (void *) &(cm->stats->hfiA[i]->always_better_than_Smax), sizeof(int),   1, fp)) goto FAILURE;
	  ESL_ALLOC(cm->stats->hfiA[i]->cm_E_cut,  sizeof(float) * cm->stats->hfiA[i]->ncut);
	  ESL_ALLOC(cm->stats->hfiA[i]->fwd_E_cut, sizeof(float) * cm->stats->hfiA[i]->ncut);
	  if (! tagged_fread(CMIO_FTHRCMECUT,  (void *)  cm->stats->hfiA[i]->cm_E_cut,                 sizeof(float), cm->stats->hfiA[i]->ncut, fp)) goto FAILURE;
	  if (! tagged_fread(CMIO_FTHRFWDECUT, (void *)  cm->stats->hfiA[i]->fwd_E_cut,                sizeof(float), cm->stats->hfiA[i]->ncut, fp)) goto FAILURE;
	  cm->stats->hfiA[i]->is_valid = TRUE; /* set valid flag */
	}
      cm->flags |= CMH_FILTER_STATS;
    }
  /* if we have exp tail stats we must have filter thresholds stats, 
   * and if we have filter threshold stats we must have exp tail stats.
   */
  if(  (cm->flags & CMH_EXPTAIL_STATS)  && (!(cm->flags & CMH_FILTER_STATS))) goto FAILURE;
  if((!(cm->flags & CMH_EXPTAIL_STATS)) &&   (cm->flags & CMH_FILTER_STATS))  goto FAILURE;

  /* Main model section */
  CMZero(cm);
  if (! tagged_fread(CMIO_STTYPE,       (void *) cm->sttype,     sizeof(char),  cm->M, fp))         goto FAILURE;
  if (! tagged_fread(CMIO_NDIDX,        (void *) cm->ndidx,      sizeof(int),   cm->M, fp))         goto FAILURE;  
  if (! tagged_fread(CMIO_STID,         (void *) cm->stid,       sizeof(char),  cm->M, fp))         goto FAILURE;  
  if (! tagged_fread(CMIO_CFIRST,       (void *) cm->cfirst,     sizeof(int),   cm->M, fp))         goto FAILURE; 
  if (! tagged_fread(CMIO_CNUM,         (void *) cm->cnum,       sizeof(int),   cm->M, fp))         goto FAILURE;
  if (! tagged_fread(CMIO_PLAST,        (void *) cm->plast,      sizeof(int),   cm->M, fp))         goto FAILURE; 
  if (! tagged_fread(CMIO_PNUM,         (void *) cm->pnum,       sizeof(int),   cm->M, fp))         goto FAILURE;  
  if (! tagged_fread(CMIO_NODEMAP,      (void *) cm->nodemap,    sizeof(int),   cm->nodes, fp))     goto FAILURE;
  if (! tagged_fread(CMIO_NDTYPE,       (void *) cm->ndtype,     sizeof(char),  cm->nodes, fp))     goto FAILURE;
  for (v = 0; v < cm->M; v++) {
    if (! tagged_fread(CMIO_T, (void *) cm->t[v], sizeof(float), MAXCONNECT, fp)) goto FAILURE;
    if (! tagged_fread(CMIO_E, (void *) cm->e[v], sizeof(float), cm->abc->K*cm->abc->K, fp)) goto FAILURE;
  }

  if (! tagged_fread(CMIO_END_DATA, (void *) NULL, 0, 0, fp)) goto FAILURE;

  /* EPN 10.29.06 Remove the sole source of CM ambiguities. Find and detach insert states
   *              that are 1 state before an END_E.  */
  cm_find_and_detach_dual_inserts(cm, 
				  FALSE, /* Don't check END_E-1 states have 0 counts, they may not if 
					  * an old version (0.7 or earlier) of cmbuild was used, or  
					  * cmbuild --nodetach  was used to build the CM  */
				  TRUE); /* Detach the states by setting trans probs into them as 0.0   */

  /* Success.
   * Renormalize the CM, and return.
   */
  CMRenormalize(cm);

  if (*ret_abc == NULL) *ret_abc = abc;	/* pass our new alphabet back to caller, if caller didn't know it already */
  *ret_cm = cm;
  return eslOK;

 FAILURE:
  if (cm != NULL) FreeCM(cm);
  ESL_FAIL(eslEFORMAT, errbuf, "Error reading the cmfile. Is it corrupt or built with a pre-1.0 cmbuild?"); 
  return status; /* NEVERREACHED */

 ERROR: 
  ESL_FAIL(eslEMEM, errbuf, "Error ran out of memory reading the cmfile.");
  return status; /* NEVERREACHED */
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
  int status;
  int tag;
  int nbytes;
  char *s;

  fread(&tag, sizeof(int), 1, fp);
  if (tag != expected_tag) return 0;
  fread(&nbytes, sizeof(int), 1, fp);
  if (nbytes > 0) {
    ESL_ALLOC(s, sizeof(char) * (nbytes+1));
    s[nbytes] = '\0';
    fread(s, sizeof(char), nbytes, fp);
  } else s = NULL; 
  *ret_s = s;
  return 1;

 ERROR:
  cm_Fail("Memory allocation error.");
  return 0; /* never reached */
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

/* EPN, Tue Aug  7 15:54:15 2007
 * is_integer() and is_real(), savagely ripped verbatim out
 * of Easel's esl_getopts.c, where they were private.
 */
/* Function: is_integer()
 * 
 * Returns TRUE if <s> points to something that atoi() will parse
 * completely and convert to an integer.
 */
static int
is_integer(char *s)
{
  int hex = 0;

  if (s == NULL) return 0;
  while (isspace((int) (*s))) s++;      /* skip whitespace */
  if (*s == '-' || *s == '+') s++;      /* skip leading sign */
				        /* skip leading conversion signals */
  if ((strncmp(s, "0x", 2) == 0 && (int) strlen(s) > 2) ||
      (strncmp(s, "0X", 2) == 0 && (int) strlen(s) > 2))
    {
      s += 2;
      hex = 1;
    }
  else if (*s == '0' && (int) strlen(s) > 1)
    s++;
				/* examine remainder for garbage chars */
  if (!hex)
    while (*s != '\0')
      {
	if (!isdigit((int) (*s))) return 0;
	s++;
      }
  else
    while (*s != '\0')
      {
	if (!isxdigit((int) (*s))) return 0;
	s++;
      }
  return 1;
}


/* is_real()
 * 
 * Returns TRUE if <s> is a string representation
 * of a valid floating point number, convertable
 * by atof().
 */
static int
is_real(char *s)
{
  int gotdecimal = 0;
  int gotexp     = 0;
  int gotreal    = 0;

  if (s == NULL) return 0;

  while (isspace((int) (*s))) s++; /* skip leading whitespace */
  if (*s == '-' || *s == '+') s++; /* skip leading sign */

  /* Examine remainder for garbage. Allowed one '.' and
   * one 'e' or 'E'; if both '.' and e/E occur, '.'
   * must be first.
   */
  while (*s != '\0')
    {
      if (isdigit((int) (*s))) 	gotreal++;
      else if (*s == '.')
	{
	  if (gotdecimal) return 0; /* can't have two */
	  if (gotexp) return 0;     /* e/E preceded . */
	  else gotdecimal++;
	}
      else if (*s == 'e' || *s == 'E')
	{
	  if (gotexp) return 0;	/* can't have two */
	  else gotexp++;
	}
      else if (isspace((int) (*s)))
	break;
      s++;
    }

  while (isspace((int) (*s))) s++;         /* skip trailing whitespace */
  if (*s == '\0' && gotreal) return 1;
  else return 0;
}

