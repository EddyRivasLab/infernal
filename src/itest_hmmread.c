/* Integration test: threaded reads of a pressed INFERNAL1/e CM database.
 *
 * A threaded cmscan reads its pressed database from several threads
 * through one CM_FILE. A thread reads a CM with cm_file_Read() while
 * holding the CM file's lock (pli_scan_mode_read_cm()), and other
 * threads read filter HMMs with cm_p7_hmmfile_Read() while holding
 * only the HMM file's lock (pli_p7_env_def()). Both use the same HMM
 * file handle, <cmfp->hfp>.
 *
 * A 1/e CM with an HMM-only HMM stores two HMMs, the filter HMM then
 * the HMM-only HMM, and cm_file_Read() reads both. It used to find the
 * second at the HMM file's current position, after the first read had
 * released the HMM file's lock. A filter HMM read by another thread in
 * between moved that position, and cm_file_Read() failed with "bad
 * file format for HMM filter" or read the wrong HMM (brief
 * 26_0821-168; seen with cmscan --cpu 8 on a whole upgraded Rfam
 * library).
 *
 * This test reads each CM once without threads to get its expected
 * filter and HMM-only HMMs (length and checksum). Then <n> threads do
 * what cmscan's threads do, in a tight loop with no search in between:
 * one reads CMs under the CM file's lock, and the others read filter
 * HMMs until it is done. Every read must succeed and give the expected
 * HMMs.
 *
 * Usage: itest_hmmread [-options] <pressed 1/e CM file>
 * The testsuite runs it on an upgraded, pressed tutorial/minifam.cm
 * (itest20-threaded-hmmread.pl).
 *
 * EPN, Sun Nov  1 2026 (w/Claude)
 */
#include <esl_config.h>
#include <p7_config.h>
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef HMMER_THREADS
#include <pthread.h>
#endif

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"

#include "infernal.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                     docgroup*/
  { "-h",        eslARG_NONE,    FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",          0 },
  { "-n",        eslARG_INT,       "4", NULL, "n>1", NULL,  NULL, NULL, "number of threads (one reads CMs, the rest HMMs)", 0 },
  { "-N",        eslARG_INT,    "1000", NULL, "n>0", NULL,  NULL, NULL, "number of CM reads",                            0 },
  { "-s",        eslARG_INT,      "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                 0 },
  { 0,0,0,0,0,0,0,0,0,0 },
};
static char usage[]  = "[-options] <pressed INFERNAL1/e CM file>";
static char banner[] = "integration test: threaded reads of a pressed CM database";

typedef struct {
  CM_FILE      *cmfp;
  ESL_ALPHABET *abc;
  int           ncm;
  off_t        *cm_offset;  /* [0..ncm-1] offset of each CM in the .i1m file           */
  off_t        *hmm_offset; /* [0..ncm-1] offset of each CM's filter HMM in the .i1m file */
  int          *fM;         /* [0..ncm-1] expected filter HMM length                    */
  uint32_t     *fck;        /* [0..ncm-1] expected filter HMM checksum                  */
  int          *hM;         /* [0..ncm-1] expected HMM-only HMM length, 0 if none       */
  uint32_t     *hck;        /* [0..ncm-1] expected HMM-only HMM checksum                */
  int           niter;      /* CM reads, for the CM reader                              */
  int           is_cm_reader;
  volatile int *done;       /* set by the CM reader when it finishes                    */
  uint32_t      seed;
  int           nfail;      /* failures seen by this thread */
  char          errmsg[2*eslERRBUFSIZE];
} WORKER_INFO;

/* hmm_checksum(): a simple checksum of an HMM's match emissions, to tell HMMs of the same length apart */
static uint32_t
hmm_checksum(const P7_HMM *hmm)
{
  uint32_t ck = 0;
  int      k, a;
  for (k = 1; k <= hmm->M; k++)
    for (a = 0; a < hmm->abc->K; a++)
      ck = ck * 31 + (uint32_t) (hmm->mat[k][a] * 1e6);
  return ck;
}

static int
read_cm_locked(WORKER_INFO *info, int i, CM_t **ret_cm)
{
  int status;
#ifdef HMMER_THREADS
  if (info->cmfp->syncRead && pthread_mutex_lock(&info->cmfp->readMutex) != 0) return eslESYS;
#endif
  cm_file_Position(info->cmfp, info->cm_offset[i]);
  status = cm_file_Read(info->cmfp, TRUE, &(info->abc), ret_cm);
#ifdef HMMER_THREADS
  if (info->cmfp->syncRead && pthread_mutex_unlock(&info->cmfp->readMutex) != 0) return eslESYS;
#endif
  return status;
}

/* worker(): the CM reader reads CMs as pli_scan_mode_read_cm() does; the
 * others read filter HMMs as pli_p7_env_def() does, until the CM
 * reader is done. (A failure message may quote another thread's
 * <cmfp->errbuf>; the count is what the test checks.)
 */
static void *
worker(void *arg)
{
  WORKER_INFO    *info = (WORKER_INFO *) arg;
  ESL_RANDOMNESS *r    = esl_randomness_Create(info->seed);
  CM_t           *cm   = NULL;
  P7_HMM         *hmm  = NULL;
  int             it, i, status;

  if (info->is_cm_reader) {
    for (it = 0; it < info->niter; it++) {
      i = esl_rnd_Roll(r, info->ncm);
      if ((status = read_cm_locked(info, i, &cm)) != eslOK) {
        if (info->nfail++ == 0) snprintf(info->errmsg, sizeof(info->errmsg), "cm_file_Read() of CM %d failed (status %d): %s", i, status, info->cmfp->errbuf);
        continue;
      }
      if (cm->fp7 == NULL || cm->fp7->M != info->fM[i] || hmm_checksum(cm->fp7) != info->fck[i] ||
          (info->hM[i] == 0 && cm->hp7 != NULL) ||
          (info->hM[i] >  0 && (cm->hp7 == NULL || cm->hp7->M != info->hM[i] || hmm_checksum(cm->hp7) != info->hck[i]))) {
        if (info->nfail++ == 0) snprintf(info->errmsg, sizeof(info->errmsg), "CM %d (%s) was read with the wrong filter or HMM-only HMM", i, cm->name);
      }
      FreeCM(cm); cm = NULL;
    }
    *(info->done) = TRUE;
  }
  else {
    while (! *(info->done)) {
      i = esl_rnd_Roll(r, info->ncm);
      if ((status = cm_p7_hmmfile_Read(info->cmfp, info->abc, info->hmm_offset[i], &hmm, NULL)) != eslOK) {
        if (info->nfail++ == 0) snprintf(info->errmsg, sizeof(info->errmsg), "cm_p7_hmmfile_Read() of HMM %d failed (status %d): %s", i, status, info->cmfp->errbuf);
        continue;
      }
      if (hmm->M != info->fM[i] || hmm_checksum(hmm) != info->fck[i]) {
        if (info->nfail++ == 0) snprintf(info->errmsg, sizeof(info->errmsg), "filter HMM %d was read wrong", i);
      }
      p7_hmm_Destroy(hmm); hmm = NULL;
    }
  }
  esl_randomness_Destroy(r);
  return NULL;
}

int
main(int argc, char **argv)
{
  ESL_GETOPTS  *go       = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char         *cmfile   = esl_opt_GetArg(go, 1);
  int           nthreads = esl_opt_GetInteger(go, "-n");
  char          errbuf[eslERRBUFSIZE];
  CM_FILE      *cmfp     = NULL;
  ESL_ALPHABET *abc      = NULL;
  P7_OPROFILE  *om       = NULL;
  CM_t         *cm       = NULL;
  WORKER_INFO  *info     = NULL;
  WORKER_INFO   shared;
  int           nalloc   = 16;
  int           nh3      = 0;
  int           nfail    = 0;
  int           i, t, status;

  if ((status = cm_file_Open(cmfile, NULL, FALSE, &cmfp, errbuf)) != eslOK) esl_fatal("Failed to open %s: %s", cmfile, errbuf);
  if (! cmfp->is_pressed)            esl_fatal("%s is not a pressed CM file (run cmpress)", cmfile);
  if (cmfp->format != CM_FILE_1e)    esl_fatal("%s is not in format INFERNAL1/e", cmfile);

  /* the offsets of each CM and its filter HMM, as cmscan gets them */
  memset(&shared, 0, sizeof(WORKER_INFO));
  ESL_ALLOC(shared.cm_offset,  sizeof(off_t)    * nalloc);
  ESL_ALLOC(shared.hmm_offset, sizeof(off_t)    * nalloc);
  while ((status = cm_p7_oprofile_ReadMSV(cmfp, TRUE, &abc, &(shared.cm_offset[shared.ncm]), NULL, NULL, NULL, NULL, NULL, &om)) == eslOK) {
    shared.hmm_offset[shared.ncm] = om->offs[p7_MOFFSET];
    p7_oprofile_Destroy(om);
    if (++shared.ncm == nalloc) {
      nalloc *= 2;
      ESL_REALLOC(shared.cm_offset,  sizeof(off_t) * nalloc);
      ESL_REALLOC(shared.hmm_offset, sizeof(off_t) * nalloc);
    }
  }
  if (status != eslEOF) esl_fatal("Failed to read the MSV profiles of %s: %s", cmfile, cmfp->msv_errbuf);
  if (shared.ncm == 0)  esl_fatal("No CMs in %s", cmfile);

  /* each CM's expected HMMs, read without threads */
  ESL_ALLOC(shared.fM,  sizeof(int)      * shared.ncm);
  ESL_ALLOC(shared.fck, sizeof(uint32_t) * shared.ncm);
  ESL_ALLOC(shared.hM,  sizeof(int)      * shared.ncm);
  ESL_ALLOC(shared.hck, sizeof(uint32_t) * shared.ncm);
  shared.cmfp = cmfp;
  shared.abc  = abc;
  for (i = 0; i < shared.ncm; i++) {
    if ((status = read_cm_locked(&shared, i, &cm)) != eslOK) esl_fatal("Failed to read CM %d: %s", i, cmfp->errbuf);
    if (cm->fp7 == NULL) esl_fatal("CM %s has no filter HMM", cm->name);
    shared.fM[i]  = cm->fp7->M;
    shared.fck[i] = hmm_checksum(cm->fp7);
    shared.hM[i]  = (cm->hp7 == NULL) ? 0 : cm->hp7->M;
    shared.hck[i] = (cm->hp7 == NULL) ? 0 : hmm_checksum(cm->hp7);
    if (cm->hp7 != NULL) nh3++;
    FreeCM(cm); cm = NULL;
  }
  if (nh3 == 0) esl_fatal("No CM in %s has an HMM-only HMM: the test needs at least one", cmfile);

#ifdef HMMER_THREADS
  {
    pthread_t *tid = NULL;
    if (cm_file_CreateLock(cmfp) != eslOK) esl_fatal("Failed to create the CM file locks");
    ESL_ALLOC(tid,  sizeof(pthread_t)   * nthreads);
    ESL_ALLOC(info, sizeof(WORKER_INFO) * nthreads);
    volatile int done = FALSE;
    for (t = 0; t < nthreads; t++) {
      info[t]       = shared;
      info[t].niter = esl_opt_GetInteger(go, "-N");
      info[t].is_cm_reader = (t == 0) ? TRUE : FALSE;
      info[t].done  = &done;
      info[t].seed  = esl_opt_GetInteger(go, "-s") + t;
      if (pthread_create(&tid[t], NULL, worker, &info[t]) != 0) esl_fatal("pthread_create() failed");
    }
    for (t = 0; t < nthreads; t++) pthread_join(tid[t], NULL);
    for (t = 0; t < nthreads; t++) {
      if (info[t].nfail > 0) printf("thread %d: %d failed reads; first: %s\n", t, info[t].nfail, info[t].errmsg);
      nfail += info[t].nfail;
    }
    free(tid);
  }
#else
  printf("itest_hmmread: not built with threads; nothing to test\n");
#endif

  free(info);
  free(shared.cm_offset); free(shared.hmm_offset);
  free(shared.fM); free(shared.fck); free(shared.hM); free(shared.hck);
  cm_file_Close(cmfp);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  if (nfail > 0) { printf("FAIL: %d failed reads\n", nfail); return 1; }
  printf("ok\n");
  return 0;

 ERROR:
  esl_fatal("allocation failed");
}
