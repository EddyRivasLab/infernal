/* cmio.c
 * SRE, Thu Aug  3 11:53:34 2000 [St. Louis]
 * CVS $Id$
 * 
 * Input/output of covariance models.
 * 
 *****************************************************************
 * @LICENSE@
 ***************************************************************** 
 */ 

#include <stdio.h>
#include "squid.h"
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
  s = MallocOrDie(sizeof(char) * (nbytes+1));
  s[nbytes] = '\0';
  fread(s, sizeof(char), nbytes, fp);
  return 1;
}

/* Function: WriteBinaryCM()
 * Date:     SRE, Thu Aug  3 12:05:30 2000 [St. Louis]
 *
 * Purpose:  Write a CM in binary format.
 *
 */
void
WriteBinaryCM(FILE *fp, CM_t *cm)
{
  int v;
  int tag;

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
  tagged_fwrite(CMIO_NDTYPE,       cm->ndtype,     sizeof(int),   cm->nodes, fp);

  for (v = 0; v < cm->M; v++) {
    tagged_fwrite(CMIO_T, cm->t[v], sizeof(float), MAXCONNECT, fp);
    tagged_fwrite(CMIO_E, cm->e[v], sizeof(float), Alphabet_size*Alphabet_size, fp);
  }

  tagged_fwrite(CMIO_END_DATA, NULL, 0, 0, fp);
}


/* Function: ReadBinaryCM()
 * Date:     SRE, Thu Aug  3 13:39:09 2000 [St. Louis]
 *
 * Purpose:  Read a CM from disk.
 */
int
ReadBinaryCM(FILE *fp, CM_t **ret_cm)
{
  CM_t         *cm;
  unsigned int  magic;
  int           M;
  int           nodes;
  int           tag;
  int           atype;
  int           asize;
  int           v;

  cm = NULL;
  if (feof(fp)) return 0;
  if (! fread((char *) &magic, sizeof(unsigned int), 1, fp)) return 0;
  if (magic != v01magic) goto FAILURE;
  
  if (! fread((char *) &M,     sizeof(int), 1, fp)) goto FAILURE;
  if (! fread((char *) &nodes, sizeof(int), 1, fp)) goto FAILURE;
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
  if (! tagged_fread(CMIO_NDTYPE,       (void *) cm->ndtype,     sizeof(int),   cm->nodes, fp))     goto FAILURE;

  for (v = 0; v < cm->M; v++) {
    if (! tagged_fread(CMIO_T, (void *) cm->t[v], sizeof(float), MAXCONNECT, fp)) goto FAILURE;
    if (! tagged_fread(CMIO_E, (void *) cm->e[v], sizeof(float), Alphabet_size*Alphabet_size, fp)) goto FAILURE;
  }
  if (! tagged_fread(CMIO_END_DATA, (void *) NULL, 0, 0, fp)) goto FAILURE;

  *ret_cm = cm;
  return 1;

 FAILURE:
  if (cm != NULL) FreeCM(cm);
  *ret_cm = NULL;
  return 1;
}
