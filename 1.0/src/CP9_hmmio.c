/* Reading/writing of CP9 HMMs no longer supported. */
#if 0
/************************************************************
 * @LICENSE@
 ************************************************************/

/* CP9_hmmio.c
 * EPN
 * 
 * Input/output of CM plan 9 HMMs.
 *
 * All following notes and functions are from HMMER 2.4's hmmio.c
 * with necessary changes for the CM plan 9 architecture. Some
 * functions from hmmio.c are omitted.
 *
 * These functions are currently never used as the ability to
 * output and input CM Plan 9 HMMs was abandoned. Currently, 
 * they're built in cmsearch or cmbuild directly from the CM
 * fast enough that I don't think it's worth inputting them
 * (it takes on the order of hundredths of seconds to build
 *  an SSU CP9 HMM). These functions did work at one point though,
 * before abandonment. Kept here for reference.
 *
 ************************************************************
 * NOTES FROM hmmio.c: 
 * As of HMMER 2.0, HMMs are saved by default in a tabular ASCII format
 * as log-odds or log probabilities scaled to an integer. A binary save
 * file format is also available which is faster to access (a
 * consideration which might be important for HMM library applications).
 * HMMs can be concatenated into HMM libraries.
 * 
 * A comment on loss of accuracy. Storing a number as a scaled log
 * probability guarantees us an error of about 0.035% or
 * less in the retrieved probability. We are relatively invulnerable
 * to the truncation errors which HMMER 1.8 was vulnerable to.  
 * 
 * Magic numbers (both for the ASCII and binary save formats) are used 
 * to label save files with a major version number. This simplifies the task of
 * backwards compatibility as new versions of the program are created. 
 * Reverse but not forward compatibility is guaranteed. I.e. HMMER 2.0
 * can read `1.7' save files, but not vice versa. Note that the major
 * version number in the save files is NOT the version of the software
 * that generated it; rather, the number of the last major version in which
 * save format changed.
 * 
 ****************************************************************** 
 * 
 * The CM Plan 9 HMM input API:
 * 
 *       CP9HMMFILE     *hmmfp;
 *       char           *hmmfile;
 *       struct cplan9_s *hmm;
 *       char            env[] = "HMMERDB";  (a la BLASTDB) 
 *
 *       hmmfp = CP9_HMMFileOpen(hmmfile, env)   NULL on failure
 *       while (CP9_HMMFileRead(hmmfp, &hmm))    0 if no more HMMs
 *          if (hmm == NULL) Die();          NULL on file parse failure
 *          whatever;
 *          FreeHMM(hmm);
 *       }
 *       CP9_HMMFileClose(hmmfp);
 *       
 *****************************************************************
 *
 * The HMM output API:
 * 
 *       FILE           *ofp;
 *       struct cplan9_s *hmm;
 *       
 *       CP9_WriteAscHMM(ofp, hmm);    to write/append an HMM to open file
 *   or  CP9_WriteBinHMM(ofp, hmm);    to write/append binary format HMM to open file
 * 
 ***************************************************************** 
 * 
 * V1.0: original implementation
 * V1.1: regularizers removed from model structure
 * V1.7: ref and cs annotation lines added from alignment, one 
 *       char per match state 1..M
 * V1.9: null model and name added to HMM structure. ASCII format changed
 *       to compact tabular one.
 * V2.0: Plan7. Essentially complete rewrite.
 */

#include "squidconf.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <unistd.h> /* to get SEEK_CUR definition on silly Suns */

#include "squid.h"
#include "funcs.h"
#include "structs.h"
#include "ssi.h"
#include "cplan9.h"

/* Magic numbers identifying binary formats.
 */
static unsigned int  vCP9magic = 0xe8ededb6; /* V2.0 binary: "hmm6" + 0x80808080 */
static unsigned int  vCP9swap  = 0xb6edede8; /* V2.0 binary, byteswapped         */

/* Old HMMER 1.x file formats.
 */
#define HMMER1_0B  1            /* binary HMMER 1.0     */
#define HMMER1_0F  2            /* flat ascii HMMER 1.0 */
#define HMMER1_1B  3            /* binary HMMER 1.1     */
#define HMMER1_1F  4            /* flat ascii HMMER 1.1 */
#define HMMER1_7B  5            /* binary HMMER 1.7     */
#define HMMER1_7F  6            /* flat ascii HMMER 1.7 */
#define HMMER1_9B  7            /* HMMER 1.9 binary     */
#define HMMER1_9F  8            /* HMMER 1.9 flat ascii */

static int  CP9_read_asc_hmm(CP9HMMFILE *hmmfp, struct cplan9_s **ret_hmm);
static int  CP9_read_bin_hmm(CP9HMMFILE *hmmfp, struct cplan9_s **ret_hmm);

static void  byteswap(char *swap, int nbytes);
static char *prob2ascii(float p, float null);
static float ascii2prob(char *s, float null);
static void  write_bin_string(FILE *fp, char *s);
static int   read_bin_string(FILE *fp, int doswap, char **ret_s);
static void  multiline(FILE *fp, char *pfx, char *s);

/*****************************************************************
 * HMM input API functions:
 *   HMMFileOpen()
 *   HMMFileRead()
 *   HMMFileClose()
 *   HMMFileRewind()
 *****************************************************************/   

/* Function: CP9_HMMFileOpen()
 * 
 * Purpose:  Open an HMM file for reading. The file may be either
 *           an index for a library of HMMs, or an HMM. 
 *           
 * Args:     hmmfile - name of file
 *           env     - NULL, or environment variable for HMM database.
 *           
 * Return:   Valid HMMFILE *, or NULL on failure.
 */
CP9HMMFILE * 
CP9_HMMFileOpen(char *hmmfile, char *env)
{
  CP9HMMFILE     *hmmfp;
  unsigned int magic;
  char         buf[512];
  char        *ssifile;
  char        *dir;        /* dir name in which HMM file was found */
  int          status;

  hmmfp = (CP9HMMFILE *) MallocOrDie (sizeof(CP9HMMFILE));
  hmmfp->f          = NULL; 
  hmmfp->parser     = NULL;
  hmmfp->is_binary  = FALSE;
  hmmfp->byteswap   = FALSE;
  hmmfp->is_seekable= TRUE;	/* always; right now, an HMM must always be in a file. */
  
  /* Open the file. Look in current directory.
   * If that doesn't work, check environment var for
   * a second possible directory (usually the location
   * of a system-wide HMM library).
   * Using dir name if necessary, construct correct SSI file name.
   */
  hmmfp->f   = NULL;
  hmmfp->ssi = NULL;
  if ((hmmfp->f = fopen(hmmfile, "r")) != NULL)
    {
      ssifile = MallocOrDie(sizeof(char) * (strlen(hmmfile) + 5));
      sprintf(ssifile, "%s.ssi", hmmfile);

      if ((hmmfp->mode = SSIRecommendMode(hmmfile)) == -1)
	Die("SSIRecommendMode() failed");
    }
  else if ((hmmfp->f = EnvFileOpen(hmmfile, env, &dir)) != NULL)
    {
      char *full;
      full    = FileConcat(dir, hmmfile);

      ssifile = MallocOrDie(sizeof(char) * (strlen(full) + strlen(hmmfile) + 5));
      sprintf(ssifile, "%s.ssi", full);

      if ((hmmfp->mode = SSIRecommendMode(full)) == -1)
	Die("SSIRecommendMode() failed");

      free(full);
      free(dir);
    }
  else return NULL;
  
  /* Open the SSI index file. If it doesn't exist, or it's corrupt, or 
   * some error happens, hmmfp->ssi stays NULL.
   */
  SQD_DPRINTF1(("Opening ssifile %s...\n", ssifile));
  SSIOpen(ssifile, &(hmmfp->ssi));
  free(ssifile);

  /* Initialize the disk offset stuff.
   */
  status = SSIGetFilePosition(hmmfp->f, hmmfp->mode, &(hmmfp->offset));
  if (status != 0) Die("SSIGetFilePosition() failed");

  /* Check for binary or byteswapped binary format
   * by peeking at first 4 bytes.
   */ 
  if (! fread((char *) &magic, sizeof(unsigned int), 1, hmmfp->f)) {
    CP9_HMMFileClose(hmmfp);
    return NULL;
  }
  rewind(hmmfp->f);

  if (magic == vCP9magic) { 
    hmmfp->parser    = CP9_read_bin_hmm;
    hmmfp->is_binary = TRUE;
    return hmmfp;
  } 
  else if (magic == vCP9swap) { 
    SQD_DPRINTF1(("Opened an Infernal CP9 HMM binary file [byteswapped]\n"));
    hmmfp->parser    = CP9_read_bin_hmm;
    hmmfp->is_binary = TRUE;
    hmmfp->byteswap  = TRUE;
    return hmmfp;
  }
  /* else we fall thru; it may be an ASCII file. */

  /* If magic looks binary but we don't recognize it, choke and die.
   */
  if (magic & 0x80000000) {
    Warn("\
%s appears to be a binary but not a CM plan 9 format that we recognize\n\
It may be from HMMER,\n\
or may be a different kind of binary altogether.\n", hmmfile);
    CP9_HMMFileClose(hmmfp);
    return NULL;
  }

  /* Check for ASCII format by peeking at first word.
   */
  if (fgets(buf, 512, hmmfp->f) == NULL) {
    CP9_HMMFileClose(hmmfp);
    return NULL;
  }
  rewind(hmmfp->f);
  
  if        (strncmp("INFERNAL-1", buf, 10) == 0) {
    hmmfp->parser = CP9_read_asc_hmm;
    return hmmfp;
  } 
  
  /* If we haven't recognized it yet, it's bogus.
   */
  CP9_HMMFileClose(hmmfp);
  return NULL;
}
int
CP9_HMMFileRead(CP9HMMFILE *hmmfp, struct cplan9_s **ret_hmm)
{
  int status;
				/* Set the disk position marker. */
  if (hmmfp->is_seekable) {
    status = SSIGetFilePosition(hmmfp->f, hmmfp->mode, &(hmmfp->offset));
    if (status != 0) Die("SSIGetFilePosition() failed");
  }
				/* Parse the HMM and return it. */
  return (*hmmfp->parser)(hmmfp, ret_hmm);
}
void
CP9_HMMFileClose(CP9HMMFILE *hmmfp)
{
  if (hmmfp->f   != NULL)  fclose(hmmfp->f);      
  if (hmmfp->ssi != NULL)  SSIClose(hmmfp->ssi);
  free(hmmfp);
}
void 
CP9_HMMFileRewind(CP9HMMFILE *hmmfp)
{
  rewind(hmmfp->f);
}
int
CP9_HMMFilePositionByName(CP9HMMFILE *hmmfp, char *name)
{	
  SSIOFFSET  offset;		/* offset in hmmfile, from SSI */
  int        fh;		/* ignored.                    */

  if (hmmfp->ssi == NULL) return 0;
  if (SSIGetOffsetByName(hmmfp->ssi, name, &fh, &offset) != 0) return 0;
  if (SSISetFilePosition(hmmfp->f, &offset) != 0) return 0;
  return 1;
}
int 
CP9_HMMFilePositionByIndex(CP9HMMFILE *hmmfp, int idx)
{				/* idx runs from 0..nhmm-1 */
  int        fh;		/* file handle is ignored; only one HMM file */
  SSIOFFSET  offset;		/* file position of HMM */

  if (hmmfp->ssi == NULL) return 0;
  if (SSIGetOffsetByNumber(hmmfp->ssi, idx, &fh, &offset) != 0) return 0;
  if (SSISetFilePosition(hmmfp->f, &offset) != 0) return 0;
  return 1;
}

/*****************************************************************
 * CP9 HMM output API:
 *    CP9_WriteAscHMM()
 *    CP9_WriteBinHMM()
 * 
 *****************************************************************/ 

/* Function: CP9_WriteAscHMM()
 * 
 * Purpose:  Save an HMM in flat text ASCII format.
 *
 * Args:     fp        - open file for writing
 *           hmm       - HMM to save
 */
void
CP9_WriteAscHMM(FILE *fp, struct cplan9_s *hmm)
{
  int k;                        /* counter for nodes             */
  int x;                        /* counter for symbols           */
  int ts;			/* counter for state transitions */

  fprintf(fp, "INFERNAL-1 [%s]\n", PACKAGE_VERSION);

  /* write header information
   */
  fprintf(fp, "NAME  %s\n", hmm->name);
  if (hmm->flags & CPLAN9_ACC)
    fprintf(fp, "ACC   %s\n", hmm->acc);
  if (hmm->flags & CPLAN9_DESC) 
    fprintf(fp, "DESC  %s\n", hmm->desc);
  fprintf(fp, "LENG  %d\n", hmm->M);
  fprintf(fp, "ALPH  %s\n", "Nucleic");   
  fprintf(fp, "RF    %s\n", (hmm->flags & CPLAN9_RF)  ? "yes" : "no");
  fprintf(fp, "CS    %s\n", (hmm->flags & CPLAN9_CS)  ? "yes" : "no");
  multiline(fp, "COM   ", hmm->comlog);
  fprintf(fp, "NSEQ  %d\n", hmm->nseq);
  fprintf(fp, "DATE  %s\n", hmm->ctime); 
  fprintf(fp, "CKSUM %d\n", hmm->checksum);
  if (hmm->flags & CPLAN9_GA)
    fprintf(fp, "GA    %.1f %.1f\n", hmm->ga1, hmm->ga2);
  if (hmm->flags & CPLAN9_TC)
    fprintf(fp, "TC    %.1f %.1f\n", hmm->tc1, hmm->tc2);
  if (hmm->flags & CPLAN9_NC)
    fprintf(fp, "NC    %.1f %.1f\n", hmm->nc1, hmm->nc2);

  /* No Specials
   */

  /* Save the null model first, so HMM readers can decode
   * log odds scores on the fly. Save as log odds probabilities
   * relative to 1/Alphabet_size (flat distribution)
   */
  fprintf(fp, "NULT  ");
  fprintf(fp, "%6s ", prob2ascii(hmm->p1, 1.0)); /* p1 */
  fprintf(fp, "%6s\n", prob2ascii(1.0-hmm->p1, 1.0));   /* p2 */
  fputs("NULE  ", fp);
  for (x = 0; x < Alphabet_size; x++)
    fprintf(fp, "%6s ", prob2ascii(hmm->null[x], 1/(float)(Alphabet_size)));
  fputs("\n", fp);

  /* EVD statistics */
  if (hmm->flags & CPLAN9_STATS) 
    fprintf(fp, "EVD   %10f %10f\n", hmm->mu, hmm->lambda);
     
  /* Print header */
  fprintf(fp, "HMM      ");
  for (x = 0; x < Alphabet_size; x++) fprintf(fp, "  %c    ", Alphabet[x]);
  fprintf(fp, "\n");
  fprintf(fp, "       %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s\n",
          "m->m", "m->i", "m->d", "i->m", "i->i", "i->d", "d->m", "d->i", "d->d", "b->m", "m->e");

  /* Print HMM parameters (main section of the save file)
   */
  for (k = 0; k <= hmm->M; k++)
    {
				/* Line 1: k, match emissions, map */
      fprintf(fp, " %5d ", k);
      for (x = 0; x < Alphabet_size; x++) 
        fprintf(fp, "%6s ", prob2ascii(hmm->mat[k][x], hmm->null[x]));
      if (hmm->flags & CPLAN9_MAP) fprintf(fp, "%5d", hmm->map[k]);
      fputs("\n", fp);
				/* Line 2: RF and insert emissions */
      fprintf(fp, " %5c ", hmm->flags & CPLAN9_RF ? hmm->rf[k] : '-');
      for (x = 0; x < Alphabet_size; x++) 
	fprintf(fp, "%6s ", (k < hmm->M) ? prob2ascii(hmm->ins[k][x], hmm->null[x]) : "*");
      fputs("\n", fp);
				/* Line 3: CS and transition probs */
      fprintf(fp, " %5c ", hmm->flags & CPLAN9_CS ? hmm->cs[k] : '-');
      for (ts = 0; ts < 9; ts++)
	fprintf(fp, "%6s ", prob2ascii(hmm->t[k][ts], 1.0)); 
      if(k > 0) fprintf(fp, "%6s ", prob2ascii(hmm->begin[k], 1.0));
      if(k > 0) fprintf(fp, "%6s ", prob2ascii(hmm->end[k], 1.0));
      
      fputs("\n", fp);
    }
  fputs("//\n", fp);
}

/* Function: CP9_WriteBinHMM()
 * 
 * Purpose:  Write a CP9 HMM in binary format.
 */
void
CP9_WriteBinHMM(FILE *fp, struct cplan9_s *hmm)
{
  int k;

  /* ye olde magic number */
  fwrite((char *) &(vCP9magic), sizeof(unsigned int), 1, fp);

  /* header section
   */
  fwrite((char *) &(hmm->flags),    sizeof(int),  1,   fp);
  write_bin_string(fp, hmm->name);
  if (hmm->flags & CPLAN9_ACC)  write_bin_string(fp, hmm->acc);
  if (hmm->flags & CPLAN9_DESC) write_bin_string(fp, hmm->desc);
  fwrite((char *) &(hmm->M),        sizeof(int),  1,   fp);
  fwrite((char *) &(Alphabet_type), sizeof(int),  1,   fp);
  if (hmm->flags & CPLAN9_RF)   fwrite((char *) hmm->rf,  sizeof(char), hmm->M+1, fp);
  if (hmm->flags & CPLAN9_CS)   fwrite((char *) hmm->cs,  sizeof(char), hmm->M+1, fp);
  if (hmm->flags & CPLAN9_MAP)  fwrite((char *) hmm->map, sizeof(int), hmm->M+1, fp);
  write_bin_string(fp, hmm->comlog);
  fwrite((char *) &(hmm->nseq),     sizeof(int),  1,   fp);
  write_bin_string(fp, hmm->ctime);
  fwrite((char *) &(hmm->checksum), sizeof(int),  1,   fp);
  if (hmm->flags & CPLAN9_GA) {
    fwrite((char *) &(hmm->ga1), sizeof(float), 1, fp);
    fwrite((char *) &(hmm->ga2), sizeof(float), 1, fp);
  }
  if (hmm->flags & CPLAN9_TC) {
    fwrite((char *) &(hmm->tc1), sizeof(float), 1, fp);
    fwrite((char *) &(hmm->tc2), sizeof(float), 1, fp);
  }
  if (hmm->flags & CPLAN9_NC) {
    fwrite((char *) &(hmm->nc1), sizeof(float), 1, fp);
    fwrite((char *) &(hmm->nc2), sizeof(float), 1, fp);
  }

  /* No Specials */

  /* Null model */
  fwrite((char *)&(hmm->p1), sizeof(float), 1,             fp);
  fwrite((char *) hmm->null, sizeof(float), Alphabet_size, fp);

  /* EVD stats */
  if (hmm->flags & CPLAN9_STATS) {
    fwrite((char *) &(hmm->mu),      sizeof(float),  1,   fp); 
    fwrite((char *) &(hmm->lambda),  sizeof(float),  1,   fp); 
  }

  /* entry/exit probabilities
   */
  fwrite((char *) hmm->begin, sizeof(float), hmm->M+1, fp);
  fwrite((char *) hmm->end,   sizeof(float), hmm->M+1, fp);

  /* main model
   */
  for (k = 0; k <= hmm->M; k++)
    fwrite((char *) hmm->mat[k], sizeof(float), Alphabet_size, fp);
  for (k = 0; k < hmm->M; k++)
    fwrite((char *) hmm->ins[k], sizeof(float), Alphabet_size, fp);
  for (k = 0; k < hmm->M; k++)
    fwrite((char *) hmm->t[k], sizeof(float), 9, fp);
}


/*****************************************************************
 *
 * Internal: HMM file parsers for CM Plan 9 HMMs.
 * 
 * CP9_read_{asc,bin}_hmm(HMMFILE *hmmfp, struct cplan9_s **ret_hmm)
 *
 * Upon return, *ret_hmm is an allocated CPlan9 HMM.
 * Return 0 if no more HMMs in the file (normal).
 * Return 1 and *ret_hmm = something if we got an HMM (normal) 
 * Return 1 if an error occurs (meaning "I tried to
 *   read something...") and *ret_hmm == NULL (meaning
 *   "...but it wasn't an HMM"). I know, this is a funny
 *   way to handle errors.
 * 
 *****************************************************************/

static int
CP9_read_asc_hmm(CP9HMMFILE *hmmfp, struct cplan9_s **ret_hmm) 
{
  struct cplan9_s *hmm;
  char  buffer[512];
  char *s;
  int   M;
  int   k, x;

  hmm = NULL;
  if (feof(hmmfp->f) || fgets(buffer, 512, hmmfp->f) == NULL) return 0;
  if (strncmp(buffer, "INFERNAL-1", 10) != 0)             goto FAILURE;

  /* Get the header information: tag/value pairs in any order,
   * ignore unknown tags, stop when "HMM" is reached (signaling
   * start of main model)
   */
  hmm = AllocCPlan9Shell();
  M = -1;
  while (fgets(buffer, 512, hmmfp->f) != NULL) {
    if      (strncmp(buffer, "NAME ", 5) == 0) CPlan9SetName(hmm, buffer+6);
    else if (strncmp(buffer, "ACC  ", 5) == 0) CPlan9SetAccession(hmm, buffer+6);
    else if (strncmp(buffer, "DESC ", 5) == 0) CPlan9SetDescription(hmm, buffer+6);
    else if (strncmp(buffer, "LENG ", 5) == 0) M = atoi(buffer+6);
    else if (strncmp(buffer, "NSEQ ", 5) == 0) hmm->nseq = atoi(buffer+6);
    else if (strncmp(buffer, "ALPH ", 5) == 0) 
      {				/* Alphabet type */
	s2upper(buffer+6);
	if (!(strncmp(buffer+6, "NUCLEIC", 7) == 0)) goto FAILURE;
      }
    else if (strncmp(buffer, "RF   ", 5) == 0) 
      {				/* Reference annotation present? */
	if (sre_toupper(*(buffer+6)) == 'Y') hmm->flags |= CPLAN9_RF;
      }
    else if (strncmp(buffer, "CS   ", 5) == 0) 
      {				/* Consensus annotation present? */
	if (sre_toupper(*(buffer+6)) == 'Y') hmm->flags |= CPLAN9_CS;
      }
    else if (strncmp(buffer, "MAP  ", 5) == 0) 
      {				/* Map annotation present? */
	if (sre_toupper(*(buffer+6)) == 'Y') hmm->flags |= CPLAN9_MAP;
      }
    else if (strncmp(buffer, "COM  ", 5) == 0) 
      {				/* Command line log */
	StringChop(buffer+6);
	if (hmm->comlog == NULL)
	  hmm->comlog = Strdup(buffer+6);
	else
	  {
	    hmm->comlog = ReallocOrDie(hmm->comlog, sizeof(char *) * 
				       (strlen(hmm->comlog) + 1 + strlen(buffer+6)));
	    strcat(hmm->comlog, "\n");
	    strcat(hmm->comlog, buffer+6);
	  }
      }
    else if (strncmp(buffer, "DATE ", 5) == 0) 
      {				/* Date file created */
	StringChop(buffer+6);
	hmm->ctime= Strdup(buffer+6); 
      }
    else if (strncmp(buffer, "GA   ", 5) == 0)
      {
	if ((s = strtok(buffer+6, " \t\n")) == NULL) goto FAILURE;
	hmm->ga1 = atof(s);
	if ((s = strtok(NULL, " \t\n")) == NULL) goto FAILURE;
	hmm->ga2 = atof(s);
	hmm->flags |= CPLAN9_GA;
      }
    else if (strncmp(buffer, "TC   ", 5) == 0)
      {
	if ((s = strtok(buffer+6, " \t\n")) == NULL) goto FAILURE;
	hmm->tc1 = atof(s);
	if ((s = strtok(NULL, " \t\n")) == NULL) goto FAILURE;
	hmm->tc2 = atof(s);
	hmm->flags |= CPLAN9_TC;
      }
    else if (strncmp(buffer, "NC   ", 5) == 0)
      {
	if ((s = strtok(buffer+6, " \t\n")) == NULL) goto FAILURE;
	hmm->nc1 = atof(s);
	if ((s = strtok(NULL, " \t\n")) == NULL) goto FAILURE;
	hmm->nc2 = atof(s);
	hmm->flags |= CPLAN9_NC;
      }
    else if (strncmp(buffer, "NULT ", 5) == 0) 
      {				/* Null model transitions */
	if ((s = strtok(buffer+6, " \t\n")) == NULL) goto FAILURE;
	hmm->p1 = ascii2prob(s, 1.);
	if ((s = strtok(NULL, " \t\n")) == NULL) goto FAILURE;
	hmm->p1 = hmm->p1 / (hmm->p1 + ascii2prob(s, 1.0));
      }
    else if (strncmp(buffer, "NULE ", 5) == 0) 
      {				/* Null model emissions */
	s = strtok(buffer+6, " \t\n");
	for (x = 0; x < Alphabet_size; x++) {
	  if (s == NULL) goto FAILURE;
	  hmm->null[x] = ascii2prob(s, 1./(float)Alphabet_size);    
	  s = strtok(NULL, " \t\n");
	}
      }
    else if (strncmp(buffer, "EVD  ", 5) == 0) 
      {				/* EVD parameters */
	hmm->flags |= CPLAN9_STATS;
	if ((s = strtok(buffer+6, " \t\n")) == NULL) goto FAILURE;
	hmm->mu = atof(s);
	if ((s = strtok(NULL, " \t\n")) == NULL) goto FAILURE;
	hmm->lambda = atof(s);
      }
    else if (strncmp(buffer, "CKSUM", 5) == 0) hmm->checksum = atoi(buffer+6);
    else if (strncmp(buffer, "HMM  ", 5) == 0) break;
  }

				/* partial check for mandatory fields */
  if (feof(hmmfp->f))                goto FAILURE;
  if (M < 1)                         goto FAILURE;
  if (hmm->name == NULL)             goto FAILURE;

  /* Main model section. Read as integer log odds, convert
   * to probabilities
   */
  AllocCPlan9Body(hmm, M);  
				/* skip an annotation line */
  if (fgets(buffer, 512, hmmfp->f) == NULL)  goto FAILURE;

				/* main model */
  for (k = 0; k <= hmm->M; k++) {
                                /* Line 1: k, match emissions, map */
    if (fgets(buffer, 512, hmmfp->f) == NULL)  goto FAILURE;
    if ((s = strtok(buffer, " \t\n")) == NULL) goto FAILURE;
    if (atoi(s) != k)                          goto FAILURE;
    for (x = 0; x < Alphabet_size; x++) {
      if ((s = strtok(NULL, " \t\n")) == NULL) goto FAILURE;
      hmm->mat[k][x] = ascii2prob(s, hmm->null[x]);
    }
    if (hmm->flags & CPLAN9_MAP) {
      if ((s = strtok(NULL, " \t\n")) == NULL) goto FAILURE;
      hmm->map[k] = atoi(s);
    }
				/* Line 2:  RF and insert emissions */
    if (fgets(buffer, 512, hmmfp->f) == NULL)  goto FAILURE;
    if ((s = strtok(buffer, " \t\n")) == NULL) goto FAILURE;
    if (hmm->flags & CPLAN9_RF) hmm->rf[k] = *s;
    for (x = 0; x < Alphabet_size; x++) {
      if ((s = strtok(NULL, " \t\n")) == NULL) goto FAILURE;
      hmm->ins[k][x] = ascii2prob(s, hmm->null[x]);
    }
    
				/* Line 3: CS and transitions */
    if (fgets(buffer, 512, hmmfp->f) == NULL)  goto FAILURE;
    if ((s = strtok(buffer, " \t\n")) == NULL) goto FAILURE;
    if (hmm->flags & CPLAN9_CS) hmm->cs[k] = *s;
    for (x = 0; x < 9; x++) {
      if ((s = strtok(NULL, " \t\n")) == NULL) goto FAILURE;
      hmm->t[k][x] = ascii2prob(s, 1.0);
    }
    if(k > 0)
      {
	if ((s = strtok(NULL, " \t\n")) == NULL) goto FAILURE;
	hmm->begin[k] = ascii2prob(s, 1.0);
	if ((s = strtok(NULL, " \t\n")) == NULL) goto FAILURE;
	hmm->end[k] = ascii2prob(s, 1.0);
      }
  } /* end loop over main model */

  /* Advance to record separator
   */
  while (fgets(buffer, 512, hmmfp->f) != NULL) 
    if (strncmp(buffer, "//", 2) == 0) break;

  /* Set flags and return
   */
  hmm->flags |= CPLAN9_HASPROB;	/* probabilities are valid */
  hmm->flags &= ~CPLAN9_HASBITS;	/* scores are not valid    */

  CPlan9Renormalize(hmm);	

  *ret_hmm = hmm;
  return 1;

FAILURE:
  if (hmm  != NULL) FreeCPlan9(hmm);
  *ret_hmm = NULL;
  return 1;
}


static int
CP9_read_bin_hmm(CP9HMMFILE *hmmfp, struct cplan9_s **ret_hmm)
{
   struct cplan9_s *hmm;
   int    k,x;
   int    type;
   unsigned int magic;

   hmm = NULL;

   /* Header section
    */
   if (feof(hmmfp->f))                                      return 0;
   if (! fread((char *) &magic, sizeof(unsigned int), 1, hmmfp->f)) return 0;

   if (hmmfp->byteswap) byteswap((char *)&magic, sizeof(unsigned int));
   if (magic != vCP9magic) goto FAILURE;
				/* allocate HMM shell for header info */
   hmm = AllocCPlan9Shell();
				/* flags */
   if (! fread((char *) &(hmm->flags), sizeof(int), 1, hmmfp->f)) goto FAILURE;
   if (hmmfp->byteswap) byteswap((char *)&(hmm->flags), sizeof(int)); 
				/* name */
   if (! read_bin_string(hmmfp->f, hmmfp->byteswap, &(hmm->name))) goto FAILURE;

				/* optional accession */
   if ((hmm->flags & CPLAN9_ACC) &&
       ! read_bin_string(hmmfp->f, hmmfp->byteswap, &(hmm->acc))) goto FAILURE;
				/* optional description */
   if ((hmm->flags & CPLAN9_DESC) &&
       ! read_bin_string(hmmfp->f, hmmfp->byteswap, &(hmm->desc))) goto FAILURE;
				/* length of model */
   if (! fread((char *) &hmm->M,  sizeof(int), 1, hmmfp->f)) goto FAILURE;
   if (hmmfp->byteswap) byteswap((char *)&(hmm->M), sizeof(int)); 
				/* alphabet type */
   if (! fread((char *) &type, sizeof(int), 1, hmmfp->f)) goto FAILURE;
   if (hmmfp->byteswap) byteswap((char *)&type, sizeof(int)); 

				/* now allocate for rest of model */
   AllocCPlan9Body(hmm, hmm->M);

				/* optional #=RF alignment annotation */
   if ((hmm->flags & CPLAN9_RF) &&
       !fread((char *) hmm->rf, sizeof(char), hmm->M+1, hmmfp->f)) goto FAILURE;
   hmm->rf[hmm->M+1] = '\0';
				/* optional #=CS alignment annotation */
   if ((hmm->flags & CPLAN9_CS) &&
       !fread((char *) hmm->cs, sizeof(char), hmm->M+1, hmmfp->f)) goto FAILURE;
   hmm->cs[hmm->M+1]  = '\0';
				/* optional alignment map annotation */
   if ((hmm->flags & CPLAN9_MAP) &&
       !fread((char *) hmm->map, sizeof(int), hmm->M+1, hmmfp->f)) goto FAILURE;
   if (hmmfp->byteswap)
     for (k = 1; k <= hmm->M; k++)
       byteswap((char*)&(hmm->map[k]), sizeof(int));
				/* command line log */
   if (!read_bin_string(hmmfp->f, hmmfp->byteswap, &(hmm->comlog)))  goto FAILURE;
				/* nseq */
   if (!fread((char *) &(hmm->nseq),sizeof(int), 1, hmmfp->f))       goto FAILURE;
   if (hmmfp->byteswap) byteswap((char *)&(hmm->nseq), sizeof(int)); 
				/* creation time */
   if (!read_bin_string(hmmfp->f, hmmfp->byteswap, &(hmm->ctime)))   goto FAILURE;
				/* checksum */
   if (!fread((char *) &(hmm->checksum),sizeof(int), 1, hmmfp->f))       goto FAILURE;
   if (hmmfp->byteswap) byteswap((char *)&(hmm->checksum), sizeof(int)); 
     
				/* Pfam gathering thresholds */
   if (hmm->flags & CPLAN9_GA) {
     if (! fread((char *) &(hmm->ga1), sizeof(float), 1, hmmfp->f)) goto FAILURE;
     if (! fread((char *) &(hmm->ga2), sizeof(float), 1, hmmfp->f)) goto FAILURE;
     if (hmmfp->byteswap) {
       byteswap((char *) &(hmm->ga1), sizeof(float));
       byteswap((char *) &(hmm->ga2), sizeof(float));
     }
   }
				/* Pfam trusted cutoffs */
   if (hmm->flags & CPLAN9_TC) {
     if (! fread((char *) &(hmm->tc1), sizeof(float), 1, hmmfp->f)) goto FAILURE;
     if (! fread((char *) &(hmm->tc2), sizeof(float), 1, hmmfp->f)) goto FAILURE;
     if (hmmfp->byteswap) {
       byteswap((char *) &(hmm->tc1), sizeof(float));
       byteswap((char *) &(hmm->tc2), sizeof(float));
     }
   }
				/* Pfam noise cutoffs */
   if (hmm->flags & CPLAN9_NC) {
     if (! fread((char *) &(hmm->nc1), sizeof(float), 1, hmmfp->f)) goto FAILURE;
     if (! fread((char *) &(hmm->nc2), sizeof(float), 1, hmmfp->f)) goto FAILURE;
     if (hmmfp->byteswap) {
       byteswap((char *) &(hmm->nc1), sizeof(float));
       byteswap((char *) &(hmm->nc2), sizeof(float));
     }
   }

   /* No specials */
   
   /* null model */
   if (!fread((char *) &(hmm->p1),sizeof(float), 1, hmmfp->f))        goto FAILURE;
   if (!fread((char *)hmm->null,sizeof(float),Alphabet_size,hmmfp->f))goto FAILURE;

  /* EVD stats */
  if (hmm->flags & CPLAN9_STATS) {
    if (! fread((char *) &(hmm->mu),     sizeof(float), 1, hmmfp->f))goto FAILURE;
    if (! fread((char *) &(hmm->lambda), sizeof(float), 1, hmmfp->f))goto FAILURE;

    if (hmmfp->byteswap) {
      byteswap((char *)&(hmm->mu),     sizeof(float));
      byteswap((char *)&(hmm->lambda), sizeof(float));
    }
  }

   /* entry/exit probabilities
    */
   if (! fread((char *) hmm->begin, sizeof(float), hmm->M+1, hmmfp->f)) goto FAILURE;
   if (! fread((char *) hmm->end,   sizeof(float), hmm->M+1, hmmfp->f)) goto FAILURE;

				/* main model */
   for (k = 0; k <= hmm->M; k++)
     if (! fread((char *) hmm->mat[k], sizeof(float), Alphabet_size, hmmfp->f)) goto FAILURE;
   for (k = 0; k <= hmm->M; k++)
     if (! fread((char *) hmm->ins[k], sizeof(float), Alphabet_size, hmmfp->f)) goto FAILURE;
   for (k = 0; k <= hmm->M; k++)
     if (! fread((char *) hmm->t[k], sizeof(float), 9, hmmfp->f)) goto FAILURE;

  /* byteswapping
   */
  if (hmmfp->byteswap) {
    for (x = 0; x < Alphabet_size; x++) 
      byteswap((char *) &(hmm->null[x]), sizeof(float));
    byteswap((char *)&(hmm->p1),   sizeof(float));

    for (k = 0; k <= hmm->M; k++) 
      { 
	for (x = 0; x < Alphabet_size; x++) 
	  byteswap((char *)&(hmm->mat[k][x]), sizeof(float));
	for (x = 0; x < Alphabet_size; x++) 
	  byteswap((char *)&(hmm->ins[k][x]), sizeof(float));
	if(k > 0)
	  {
	    byteswap((char *)&(hmm->begin[k]),  sizeof(float));
	    byteswap((char *)&(hmm->end[k]),    sizeof(float));
	  }
	if (k < hmm->M)
	  for (x = 0; x < 9; x++) 
	    byteswap((char *)&(hmm->t[k][x]), sizeof(float));
      }
  }

    
  /* set flags and return
   */
  hmm->flags |= CPLAN9_HASPROB;	        /* probabilities are valid  */
  hmm->flags &= ~CPLAN9_HASBITS;	/* scores are not yet valid */
  *ret_hmm = hmm;
  return 1;

FAILURE:
  if (hmm != NULL) FreeCPlan9(hmm);
  *ret_hmm = NULL;
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
  static char buffer[8];

  if (p == 0.0) return "*";
  sprintf(buffer, "%6d", Prob2Score(p, null));
  return buffer;
}


/* Function: ascii2prob()
 * 
 * Purpose:  Convert a saved string back to a probability.
 */
static float
ascii2prob(char *s, float null)
{
  return (*s == '*') ? 0. : Score2Prob(atoi(s), null);
}

/* Function: byteswap()
 * 
 * Purpose:  Swap between big-endian and little-endian.
 *           For example:
 *               int foo = 0x12345678;
 *               byteswap((char *) &foo, sizeof(int));
 *               printf("%x\n", foo)
 *           gives 78563412.
 *           
 *           I don't fully understand byte-swapping issues.
 *           However, I have tested this on chars through floats,
 *           on various machines:
 *               SGI IRIX 4.0.5, SunOS 4.1.3, DEC Alpha OSF/1, Alliant
 *               
 *           Note: this is only a partial solution to the problem of
 *           binary file portability. 32 bit integers are assumed by HMMER,
 *           for instance. This should be true for all UNIX, VAX, and WinNT
 *           platforms, I believe.     
 *
 * Date: Sun Feb 12 10:26:22 1995              
 */
static void
byteswap(char *swap, int nbytes)
{
  int  x;
  char byte;
  
  for (x = 0; x < nbytes / 2; x++)
    {
      byte = swap[nbytes - x - 1];
      swap[nbytes - x - 1] = swap[x];
      swap[x] = byte;
    }
}

/* Function: write_bin_string()
 * Date:     SRE, Wed Oct 29 13:49:27 1997 [TWA 721 over Canada]
 * 
 * Purpose:  Write a string in binary save format: an integer
 *           for the string length (including \0), followed by
 *           the string.
 */
static void
write_bin_string(FILE *fp, char *s)
{
  int len;
  if (s != NULL) 
    {
      len = strlen(s) + 1;
      fwrite((char *) &len, sizeof(int),  1,   fp);
      fwrite((char *) s,    sizeof(char), len, fp);
    }
  else
    {
      len = 0;
      fwrite((char *) &len, sizeof(int), 1, fp);
    }
}

/* Function: read_bin_string()
 * Date:     SRE, Wed Oct 29 14:03:23 1997 [TWA 721]
 * 
 * Purpose:  Read in a string from a binary file, where
 *           the first integer is the length (including '\0').
 *           
 * Args:     fp       - FILE to read from
 *           doswap   - TRUE to byteswap
 *           ret_s    - string to read into
 *                             
 * Return:   0 on failure. ret_s is malloc'ed here.
 */                            
static int
read_bin_string(FILE *fp, int doswap, char **ret_s)
{
  char *s;
  int   len;

  if (! fread((char *) &len, sizeof(int), 1, fp))  return 0;
  if (doswap) byteswap((char *)&len, sizeof(int)); 
  s = MallocOrDie (sizeof(char) * (len));
  if (! fread((char *) s, sizeof(char), len, fp)) 
    {
      free(s);
      return 0;
    }

  *ret_s = s;
  return 1;
}

/* Function: multiline()
 * Date:     Mon Jan  5 14:57:50 1998 [StL]
 * 
 * Purpose:  Given a record (like the comlog) that contains 
 *           multiple lines, print it as multiple lines with
 *           a given prefix. e.g.:
 *           
 *           given:   "COM   ", "foo\nbar\nbaz"
 *           print:   COM   foo
 *                    COM   bar
 *                    COM   baz
 *                    
 *                    
 *           Used to print the command log to ASCII save files.
 *           
 * Args:     fp:   FILE to print to
 *           pfx:  prefix for each line
 *           s:    line to break up and print; tolerates a NULL
 *
 * Return:   (void)
 */
static void
multiline(FILE *fp, char *pfx, char *s)
{
  char *buf;
  char *sptr;

  if (s == NULL) return;
  buf  = Strdup(s);
  sptr = strtok(buf, "\n");
  while (sptr != NULL)
    {
      fprintf(fp, "%s%s\n", pfx, sptr);
      sptr = strtok(NULL, "\n");
    }
  free(buf);
}
#endif
