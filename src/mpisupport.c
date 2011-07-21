/* mpisupport.c
 * Optional support for MPI parallelization in Infernal.
 * 
 * EPN, Mon Aug 27 12:38:13 2007
 * SVN $Id$
 */

#include "esl_config.h"
#include "p7_config.h"
#include "config.h"

#ifdef HAVE_MPI

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "mpi.h"

#include "easel.h"
#include "esl_mpi.h"

#include "funcs.h"
#include "structs.h"

static int cm_hit_MPISend(CM_HIT *hit, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc);
static int cm_hit_MPIPackSize(CM_HIT *hit, MPI_Comm comm, int *ret_n);
static int cm_hit_MPIPack(CM_HIT *hit, char *buf, int n, int *pos, MPI_Comm comm);
static int cm_hit_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, CM_HIT *hit);
static int cm_hit_MPIRecv(int source, int tag, MPI_Comm comm, char **buf, int *nalloc, CM_HIT *hit);

/*****************************************************************
 * 1. Communicating a CM.
 *****************************************************************/

/* Function:  cm_master_MPIBcast()
 * Incept:    EPN, Wed May  9 17:24:53 2007
 *
 * Purpose:   Broadcasts CM <cm> from a master.
 *            
 *            If <cm> is NULL, broadcasts a end-of-data signal, to
 *            tell workers to shut down.
 *            
 */
int
cm_master_MPIBcast(CM_t *cm, int tag, MPI_Comm comm, char **buf, int *nalloc)
{
  int   status;
  int   code;
  int   sz, n, pos;

  /* Figure out size */
  if (MPI_Pack_size(1, MPI_INT, comm, &n) != 0) ESL_XEXCEPTION(eslESYS, "mpi pack size failed");
  if (cm != NULL) {
    if ((status = cm_justread_MPIPackSize(cm, comm, &sz)) != eslOK) return status; 
    n += sz;
  }

  /* Make sure the buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) {
    void *tmp;
    ESL_RALLOC(*buf, tmp, sizeof(char) * n);
    *nalloc = n; 
  }

  /* Pack the status code and CM into the buffer */
  pos  = 0;
  code = (cm == NULL) ? eslEOD : eslOK;
  if (MPI_Pack(&code, 1, MPI_INT, *buf, n, &pos, comm) != 0) ESL_EXCEPTION(eslESYS, "mpi pack failed");
  if (cm != NULL) {
    if ((status = cm_justread_MPIPack(cm, *buf, n, &pos, comm)) != eslOK) return status;
  }

  /* Broadcast the size of the packed CM */
  if (MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD) != 0) ESL_EXCEPTION(eslESYS, "mpi broadcast failed");
  /* Broadcast the packed CM */
  if (MPI_Bcast (*buf, n, MPI_PACKED, 0, MPI_COMM_WORLD) != 0) ESL_EXCEPTION(eslESYS, "mpi broadcast failed");
  return eslOK;
  
 ERROR:
  return status;
}

/* Function:  cm_worker_MPIBcast()
 * Incept:    EPN, Mon Aug 27 14:44:11 2007
 *
 * Purpose:   Receives a broadcasted buffer size and then a
 *            broadcasted work unit of that size that consists
 *            of a single CM <cm> tagged as <tag> from 
 *            communicator <comm> (usually 0, a master).
 *            
 *            Work units are prefixed by a status code. If the unit's
 *            code is <eslOK> and no errors are encountered, this
 *            routine will return <eslOK> and a non-<NULL> <*ret_cm>.
 *            If the unit's code is <eslEOD> (a shutdown signal), 
 *            this routine returns <eslEOD> and <*ret_cm> is <NULL>.
 *            
 *            If the packed CM is an end-of-data signal, return
 *            <eslEOD>, and <*ret_cm> is <NULL>.
 *            
 * Returns:   <eslOK> on success. <*ret_cm> contains the new CM; it
 *            is allocated here, and the caller is responsible for
 *            free'ing it.  <*buf> may have been reallocated to a
 *            larger size, and <*nalloc> may have been increased.
 *
 *
 * Throws:    <eslESYS> if an MPI call fails; <eslEMEM> if an allocation fails.
 *            In either case, <*ret_msa> is NULL, and the <buf> and its size
 *            <*nalloc> remain valid.
 * Xref:      J1/72.
 *            
 */
int
cm_worker_MPIBcast(int tag, MPI_Comm comm, char **buf, int *nalloc, ESL_ALPHABET **abc, CM_t **ret_cm)
{
  int         status, code;
  int         n;
  int         pos;

  /* Receive the buffer size broadcast */
  if (MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD) != 0) ESL_EXCEPTION(eslESYS, "mpi broadcast failed.");

  /* Make sure the buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) {
    void *tmp;
    ESL_RALLOC(*buf, tmp, sizeof(char) * n);
    *nalloc = n; 
  }

  /* Receive the CM broadcast */
  if (MPI_Bcast (*buf, n, MPI_PACKED, 0, MPI_COMM_WORLD) != 0) ESL_EXCEPTION(eslESYS, "mpi broadcast failed.");

  /* Unpack it - where the first integer is a status code, OK or EOD */
  pos = 0;
  if (MPI_Unpack       (*buf, n, &pos, &code, 1, MPI_INT, comm) != 0)  ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (code == eslEOD) { status = eslEOD; goto ERROR; }

  return cm_justread_MPIUnpack(abc, *buf, *nalloc, &pos, comm, ret_cm);

 ERROR:
  *ret_cm = NULL;
  return status;
}

/* Function:  cm_justread_MPIUnpack()
 * Synopsis:  Unpacks a CM <cm> (packed by cm_justread_MPIPack()) from an MPI buffer.
 * Incept:    EPN, Mon Aug 27 14:44:11 2007
 *
 * Purpose:   Unpack a newly allocated CM from MPI packed buffer
 *            <buf>, starting from position <*pos>, where the total length
 *            of the buffer in bytes is <n>. Note that 
 *            the full CM data structure was not packed, only
 *            the parts that could have been changed from their initial
 *            (default) values in CMFileRead() are packed. In this
 *            case it's safe to 'fill-in' the 'rest' of the CM with a ConfigCM()
 *            call, which saves us from broadcasting the 'rest' of the CM.
 * 
 * Returns:   <eslOK> on success. <*pos> is updated to the position of
 *            the next element in <buf> to unpack (if any). <*ret_cm>
 *            contains a newly allocated CM, which the caller is 
 *            responsible for free'ing.
 *            
 * Throws:    <eslESYS> on an MPI call failure. <eslEMEM> on allocation failure.
 *            In either case, <*ret_cm> is <NULL>, and the state of <buf>
 *            and <*pos> is undefined and should be considered to be corrupted.
 */
int
cm_justread_MPIUnpack(ESL_ALPHABET **abc, char *buf, int n, int *pos, MPI_Comm comm, CM_t **ret_cm)
{
  int    status;
  CM_t  *cm                = NULL;
  P7_HMM *ahmm             = NULL;
  ESL_ALPHABET *aabc       = NULL;
  int    atype, i, K, M, nnodes;
  float *tmp_ap7_gfmuA     = NULL;
  float *tmp_ap7_gflambdaA = NULL;
  int    nap7_expected = 0;
  char   errbuf[cmERRBUFSIZE];

  cm = CreateCMShell(); 
  if (cm == NULL) { status = eslEMEM; goto ERROR; }
  if (MPI_Unpack(buf, n, pos, &(cm->M),     1, MPI_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->nodes), 1, MPI_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->clen),  1, MPI_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &atype,       1, MPI_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->flags), 1,   MPI_INT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  /* note: it's important cm->flags is set prior to CreateCMBody() call below, to properly allocate
   * cm->rf, cm->consensus, cm->map, or not */

  /* Set or verify the alphabet */
  if (*abc == NULL)	{	/* still unknown: set it, pass control of it back to caller */
    if ((*abc = esl_alphabet_Create(atype)) == NULL)       { status = eslEMEM;      goto ERROR; }
  } else {			/* already known: check it */
    if ((*abc)->type != atype)                             { status = eslEINCOMPAT; goto ERROR; }
  }

  CreateCMBody(cm, cm->nodes, cm->M, cm->clen, (*abc));

  /* For convenience below */
  K      = cm->abc->K;
  M      = cm->M;
  nnodes = cm->nodes; 

  /* Unpack the rest of the CM */
  /* note cm->flags is how it was immediately after CM is read from file, so the only flags that can possibly be raised are:
   * CMH_ACC, CMH_DESC, CM_RF, CMH_GA, CMH_TC, CMH_NC, CMH_CHKSUM, CMH_MAP, CMH_CONS, CMH_EXPTAIL_STATS, CMH_FILTER_STATS, 
   * CMH_QDB, CMH_MLP7_STATS, CMH_AP7_STATS, CM_IS_RSEARCH, CM_RSEARCHTRANS, CM_RSEARCH_EMIT */

  if (MPI_Unpack(buf, n, pos, cm->e[0],              M*K*K, MPI_FLOAT, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, cm->t[0],       M*MAXCONNECT, MPI_FLOAT, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, cm->sttype,              M+1,  MPI_CHAR, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, cm->stid,                M+1,  MPI_CHAR, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, cm->ndidx,                 M,   MPI_INT, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, cm->cfirst,                M,   MPI_INT, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, cm->cnum,                  M,   MPI_INT, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, cm->plast,                 M,   MPI_INT, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, cm->pnum,                  M,   MPI_INT, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, cm->nodemap,          nnodes,   MPI_INT, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, cm->ndtype,           nnodes,  MPI_CHAR, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  if (MPI_Unpack(buf, n, pos, &(cm->nseq),               1,   MPI_INT, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->eff_nseq),           1, MPI_FLOAT, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->el_selfsc),          1, MPI_FLOAT, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->null2_omega),        1, MPI_FLOAT, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->null3_omega),        1, MPI_FLOAT, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->checksum),           1,   MPI_INT, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, cm->null,                  K, MPI_FLOAT, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  if (MPI_Unpack(buf, n, pos, cm->dmin,                  M,   MPI_INT, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, cm->dmax,                  M,   MPI_INT, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->W),                  1,   MPI_INT, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->beta_W),             1,MPI_DOUBLE, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->beta_qdb),           1,MPI_DOUBLE, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->tau),                1,MPI_DOUBLE, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->ga),                 1, MPI_FLOAT, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->nc),                 1, MPI_FLOAT, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->tc),                 1, MPI_FLOAT, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  /* free the cm->comlog that was allocated inside the CreateCMBody() call, we'll allocate a new one, a bit messy */
  FreeComLog(cm->comlog);
  if((status = comlog_MPIUnpack(buf, n, pos, comm, &(cm->comlog))) != eslOK) ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 

  /* name and all the optional stuff */
  if ((status = esl_mpi_UnpackOpt(buf, n, pos, (void**)&(cm->name), NULL, MPI_CHAR, comm)) != eslOK) ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
  if ((status = esl_mpi_UnpackOpt(buf, n, pos, (void**)&(cm->acc),  NULL, MPI_CHAR, comm)) != eslOK) ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
  if ((status = esl_mpi_UnpackOpt(buf, n, pos, (void**)&(cm->desc), NULL, MPI_CHAR, comm)) != eslOK) ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 

  if (cm->flags & CMH_RF)   { if (MPI_Unpack(buf, n, pos, cm->rf,        cm->clen+2, MPI_CHAR, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); }
  if (cm->flags & CMH_CONS) { if (MPI_Unpack(buf, n, pos, cm->consensus, cm->clen+2, MPI_CHAR, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); }
  if (cm->flags & CMH_MAP)  { if (MPI_Unpack(buf, n, pos, cm->map,       cm->clen+1, MPI_INT,  comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); }

  /* the filter stats */
  if (cm->flags & CMH_FILTER_STATS) { 
    ESL_ALLOC(cm->hfiA, sizeof(HMMFilterInfo_t *) * FTHR_NMODES);
    for(i = 0; i < FTHR_NMODES; i++) { 
      if((status = hmm_filter_info_MPIUnpack(buf, n, pos, comm, &(cm->hfiA[i]))) != eslOK) ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
    }
  }  

  /* the E-value stats */
  if (cm->flags & CMH_EXPTAIL_STATS) { 
    ESL_ALLOC(cm->expA, sizeof(ExpInfo_t *) * EXP_NMODES);
    for(i = 0; i < EXP_NMODES; i++) { 
      if((status = exp_info_MPIUnpack(buf, n, pos, comm, &(cm->expA[i]))) != eslOK) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
    }
  }
  if (cm->flags & CMH_MLP7_STATS) { 
    if (MPI_Unpack(buf, n, pos, &(cm->mlp7_evparam),  CM_p7_NEVPARAM, MPI_FLOAT, comm)  != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  }

  /* Finally unpack any additional p7 filter stats, and the p7 models themselves */
  if (MPI_Unpack(buf, n, pos, &nap7_expected,  1,   MPI_INT, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (nap7_expected > 0) { 
    if (cm->flags & CMH_AP7_STATS) { 
      ESL_ALLOC(tmp_ap7_gfmuA,     sizeof(float) * nap7_expected);
      ESL_ALLOC(tmp_ap7_gflambdaA, sizeof(float) * nap7_expected);
      if (MPI_Unpack(buf, n, pos, tmp_ap7_gfmuA,     nap7_expected, MPI_FLOAT, comm)  != 0)      ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
      if (MPI_Unpack(buf, n, pos, tmp_ap7_gflambdaA, nap7_expected, MPI_FLOAT, comm)  != 0)      ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
    }
    for(i = 0; i < nap7_expected; i++) { 
      if((status = p7_hmm_MPIUnpack(buf, n, pos, comm, &aabc, &ahmm)) != eslOK) ESL_XEXCEPTION(status, "mpi unpack failed");
      if((status = cm_Addp7(cm, ahmm, tmp_ap7_gfmuA[i], tmp_ap7_gflambdaA[i], errbuf)) != eslOK) ESL_XEXCEPTION(status, "mpi unpack failed:\n%s", errbuf);
    }      
  }

  if(tmp_ap7_gfmuA     != NULL) free(tmp_ap7_gfmuA);
  if(tmp_ap7_gflambdaA != NULL) free(tmp_ap7_gflambdaA);

  *ret_cm = cm;
  return eslOK;
  
  ERROR:
  if(cm   != NULL) FreeCM(cm);
  if(tmp_ap7_gfmuA     != NULL) free(tmp_ap7_gfmuA);
  if(tmp_ap7_gflambdaA != NULL) free(tmp_ap7_gflambdaA);
  return status;
}

/* Function:  cm_justread_MPIPack()
 * Incept:    EPN, Mon Aug 27 14:24:57 2007
 *
 * Purpose:   Packs CM <cm> into an MPI packed message buffer <buf>
 *            of length <n> bytes, starting at byte position <*position>,
 *            for MPI communicator <comm>. Note that 
 *            the full CM data structure is not packed, only
 *            the parts that could have been changed from their initial
 *            (default) values in CMFileRead() are packed. In this
 *            case it's safe to 'fill-in' the 'rest' of the CM with a ConfigCM()
 *            call, which saves us from broadcasting the 'rest' of the CM.
 *            
 *            The caller must know that <buf>'s allocation of <n>
 *            bytes is large enough to append the packed CM at
 *            position <*pos>. This typically requires a call to
 *            <cm_justread_MPIPackSize()> first, and reallocation if
 *            needed.
 *            
 * Returns:   <eslOK> on success; <buf> now contains the
 *            packed <cm>, and <*position> is set to the byte
 *            immediately following the last byte of the CM
 *            in <buf>. 
 *
 * Throws:    <eslESYS> if an MPI call fails; or <eslEMEM> if the
 *            buffer's length <n> was overflowed in trying to pack
 *            <cm> into <buf>. In either case, the state of
 *            <buf> and <*position> is undefined, and both should
 *            be considered to be corrupted.
 */
int
cm_justread_MPIPack(CM_t *cm, char *buf, int n, int *pos, MPI_Comm comm)
{
  int   status;
  int   i;
  int   K      = cm->abc->K;
  int   M      = cm->M;
  int   nnodes = cm->nodes;
  int   clen   = cm->clen;
  float *tmp_ap7_gfmuA     = NULL;
  float *tmp_ap7_gflambdaA = NULL;

  if (MPI_Pack(&M,                        1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&nnodes,                   1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&clen,                     1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack((void *) &(cm->abc->type), 1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->flags),              1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "pack failed");

  /* note cm->flags is how it was immediately after CM is read from file, so the only flags that can possibly be raised are:
   * CMH_ACC, CMH_DESC, CM_RF, CMH_GA, CMH_TC, CMH_NC, CMH_CHKSUM, CMH_MAP, CMH_CONS, CMH_EXPTAIL_STATS, CMH_FILTER_STATS, 
   * CMH_QDB, CMH_MLP7_STATS, CMH_AP7_STATS, CM_IS_RSEARCH, CM_RSEARCHTRANS, CM_RSEARCH_EMIT */

  if (MPI_Pack(cm->e[0],              M*K*K, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(cm->t[0],       M*MAXCONNECT, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(cm->sttype,              M+1,  MPI_CHAR, buf, n, pos, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(cm->stid,                M+1,  MPI_CHAR, buf, n, pos, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(cm->ndidx,                 M,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(cm->cfirst,                M,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(cm->cnum,                  M,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(cm->plast,                 M,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(cm->pnum,                  M,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(cm->nodemap,          nnodes,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(cm->ndtype,           nnodes,  MPI_CHAR, buf, n, pos, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "pack failed");

  if (MPI_Pack(&(cm->nseq),               1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->eff_nseq),           1, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->el_selfsc),          1, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->null2_omega),        1, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->null3_omega),        1, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->checksum),           1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(cm->null,                  K, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "pack failed");

  if (MPI_Pack(cm->dmin,                  M,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(cm->dmax,                  M,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->W),                  1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->beta_W),             1,MPI_DOUBLE, buf, n, pos, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->beta_qdb),           1,MPI_DOUBLE, buf, n, pos, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->tau),                1,MPI_DOUBLE, buf, n, pos, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->ga),                 1, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->nc),                 1, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->tc),                 1, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "pack failed");

  /* comlog */
  if((status = comlog_MPIPack(cm->comlog, buf, n, pos, comm)) != eslOK)               ESL_XEXCEPTION(status, "pack failed");

  /* name and all the optional stuff */
  if ((status = esl_mpi_PackOpt(cm->name,    -1, MPI_CHAR, buf, n, pos, comm)) != eslOK) ESL_XEXCEPTION(eslESYS, "pack failed");
  if ((status = esl_mpi_PackOpt(cm->acc,     -1, MPI_CHAR, buf, n, pos, comm)) != eslOK) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if ((status = esl_mpi_PackOpt(cm->desc,    -1, MPI_CHAR, buf, n, pos, comm)) != eslOK) ESL_XEXCEPTION(eslESYS, "pack failed"); 

  if (cm->flags & CMH_RF)   { if (MPI_Pack(cm->rf,        cm->clen+2, MPI_CHAR, buf, n, pos, comm)  != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); }
  if (cm->flags & CMH_CONS) { if (MPI_Pack(cm->consensus, cm->clen+2, MPI_CHAR, buf, n, pos, comm)  != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); }
  if (cm->flags & CMH_MAP)  { if (MPI_Pack(cm->map,       cm->clen+1, MPI_INT,  buf, n, pos, comm)  != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); }

  /* the filter stats */
  if (cm->flags & CMH_FILTER_STATS) { 
    for(i = 0; i < FTHR_NMODES; i++) { 
      if ((status = hmm_filter_info_MPIPack(cm->hfiA[i], buf, n, pos, comm)) != eslOK)  ESL_XEXCEPTION(eslESYS, "pack failed");
    }
  }  

  /* the E-value stats */
  if (cm->flags & CMH_EXPTAIL_STATS) { 
    for(i = 0; i < EXP_NMODES; i++) { 
      if ((status = exp_info_MPIPack(cm->expA[i], buf, n, pos, comm)) != eslOK) ESL_XEXCEPTION(eslESYS, "pack failed");
    }
  }
  if (cm->flags & CMH_MLP7_STATS) { 
    if (MPI_Pack(&(cm->mlp7_evparam), CM_p7_NEVPARAM, MPI_FLOAT, buf, n, pos, comm)  != 0) ESL_XEXCEPTION(eslESYS, "pack failed");
  }

  /* Finally pack any additional p7 filter stats, and the p7 models themselves */
  if (MPI_Pack(&(cm->nap7),               1, MPI_INT, buf, n, pos, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "pack failed");
  if(cm->nap7 > 0) { 
    if (cm->flags & CMH_AP7_STATS) { 
      ESL_ALLOC(tmp_ap7_gfmuA,     sizeof(float) * cm->nap7);
      ESL_ALLOC(tmp_ap7_gflambdaA, sizeof(float) * cm->nap7);
      for(i = 0; i < cm->nap7; i++) { 
	tmp_ap7_gfmuA[i]     = cm->ap7_evparamAA[i][CM_p7_GFMU];
	tmp_ap7_gflambdaA[i] = cm->ap7_evparamAA[i][CM_p7_GFLAMBDA];
      }
      if (MPI_Pack(tmp_ap7_gfmuA,     cm->nap7, MPI_FLOAT, buf, n, pos, comm)  != 0) ESL_XEXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack(tmp_ap7_gflambdaA, cm->nap7, MPI_FLOAT, buf, n, pos, comm)  != 0) ESL_XEXCEPTION(eslESYS, "pack failed");
    }
    for(i = 0; i < cm->nap7; i++) { 
      if((status = p7_hmm_MPIPack(cm->ap7A[i], buf, n, pos, comm)) != eslOK) ESL_XEXCEPTION(status, "pack failed");
    }
  }

  if (*pos > n) ESL_XEXCEPTION(eslEMEM, "buffer overflow");

  if(tmp_ap7_gfmuA     != NULL) free(tmp_ap7_gfmuA);
  if(tmp_ap7_gflambdaA != NULL) free(tmp_ap7_gflambdaA);

  return eslOK;

  ERROR: 
  if(tmp_ap7_gfmuA     != NULL) free(tmp_ap7_gfmuA);
  if(tmp_ap7_gflambdaA != NULL) free(tmp_ap7_gflambdaA);
  return status;
}

/* Function:  cm_justread_MPIPackSize()
 *
 * Synopsis:  Calculates size needed to pack a CM that has
 *            just been read from a CM file by a CMFileRead()
 *            call, we'll need to pack far less than a fully
 *            configure CM in this case.
 *
 * Incept:    EPN, Mon Aug 27 10:34:15 2007
 *            based on p7_hmm_MPIPackSize() from HMMER3.
 *
 * Purpose:   Calculate an upper bound on the number of bytes
 *            that <cm_justread_MPIPack()> will need to pack a CM
 *            <cm> in a packed MPI message for MPI communicator
 *            <comm>; return that number of bytes in <*ret_n>.
 *
 * Returns:   <eslOK> on success, and <*ret_n> contains the answer.
 *
 * Throws:    <eslESYS> if an MPI call fails, and <*ret_n> is 0.
 */
int
cm_justread_MPIPackSize(CM_t *cm, MPI_Comm comm, int *ret_n)
{
  int   status;
  int   i;
  int   n = 0;
  int   K = cm->abc->K;
  int   M = cm->M;
  int   nnodes = cm->nodes;
  int   sz;

  if (MPI_Pack_size(1,            MPI_INT,    comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  
  n += 5*sz; /* M, nodes, clen, abc->type, flags */ 
  if (MPI_Pack_size(M*K*K,        MPI_FLOAT,  comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  
  n += sz; /* e */
  if (MPI_Pack_size(M*MAXCONNECT, MPI_FLOAT,  comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  
  n += sz; /* t */
  if (MPI_Pack_size(M+1,          MPI_CHAR,   comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  
  n += 2*sz; /* sttype, stid */
  if (MPI_Pack_size(M,            MPI_INT,    comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  
  n += 5*sz; /* ndidx, cfirst, cnum, plast, pnum */
  if (MPI_Pack_size(nnodes,       MPI_INT,    comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  
  n += sz; /* nodemap */
  if (MPI_Pack_size(nnodes,       MPI_CHAR,   comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  
  n += sz; /* ndtype */
  if (MPI_Pack_size(1,            MPI_INT,    comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  
  n += sz; /* nseq */
  if (MPI_Pack_size(1,            MPI_FLOAT,  comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  
  n += 4*sz; /* eff_nseq, el_selfsc, null2_omega, null3_omega */
  if (MPI_Pack_size(1,            MPI_INT,    comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  
  n += sz; /* checksum */
  if (MPI_Pack_size(K,            MPI_FLOAT,  comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += sz; /* null */
  if (MPI_Pack_size(M,            MPI_INT,    comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  
  n += 2*sz; /* dmin/dmax */
  if (MPI_Pack_size(1,            MPI_INT,    comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  
  n += sz; /* W */
  if (MPI_Pack_size(1,            MPI_DOUBLE, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += 3*sz; /* beta_W, beta_qdb, tau */
  if (MPI_Pack_size(1,            MPI_FLOAT,  comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += 3*sz; /* ga, tc, nc */

  if ((status = comlog_MPIPackSize(cm->comlog, comm, &sz)) != eslOK) ESL_XEXCEPTION(eslESYS, "pack size failed");
  printf("comlog size: %d\n", sz);
  n += sz; /* comlog */

  if ((status = esl_mpi_PackOptSize(cm->name, -1, MPI_CHAR, comm, &sz)) != eslOK) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += sz; /* name */
  if ((status = esl_mpi_PackOptSize(cm->acc,  -1, MPI_CHAR, comm, &sz)) != eslOK) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += sz; /* acc */
  if ((status = esl_mpi_PackOptSize(cm->desc, -1, MPI_CHAR, comm, &sz)) != eslOK) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += sz; /* desc */

  if (cm->flags & CMH_RF)   if (MPI_Pack_size(cm->clen+2,  MPI_CHAR, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += sz; /* rf */
  if (cm->flags & CMH_CONS) if (MPI_Pack_size(cm->clen+2,  MPI_CHAR, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += sz; /* consensus */
  if (cm->flags & CMH_MAP)  if (MPI_Pack_size(cm->clen+1,   MPI_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += sz; /* map */

  if (cm->flags & CMH_FILTER_STATS) { 
    for(i = 0; i < FTHR_NMODES; i++) { 
      if ((status = hmm_filter_info_MPIPackSize(cm->hfiA[i], comm, &sz)) != eslOK) ESL_XEXCEPTION(eslESYS, "pack size failed");
      n += sz;
    }
  }

  if (cm->flags & CMH_EXPTAIL_STATS) { 
    for(i = 0; i < EXP_NMODES; i++) { 
      if ((status = exp_info_MPIPackSize(cm->expA[i], comm, &sz)) != eslOK) ESL_XEXCEPTION(eslESYS, "pack size failed");
      n += sz;
    }
  }
  if (cm->flags & CMH_MLP7_STATS) { 
    if (MPI_Pack_size(CM_p7_NEVPARAM,  MPI_FLOAT,  comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
    n += sz;
  }

  /* Finally, compute size for additional p7 filters */
  if (MPI_Pack_size(1, MPI_INT,  comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += sz; /* cm->nap7 */
  if(cm->nap7 > 0) { 
    if (MPI_Pack_size(1, MPI_FLOAT,  comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
    n += 2*cm->nap7*sz; /* tmp_ap7_gfmuA and tmp_ap7_gflambdaA */
    for(i = 0; i < cm->nap7; i++) { 
      if ((status = p7_hmm_MPIPackSize(cm->ap7A[i], comm, &sz)) != eslOK) ESL_XEXCEPTION(status, "p7_hmm_MPIPackSize() failed");
      n += sz; /* cm->ap7A[i] p7 HMM */
    }
  }

  *ret_n = n;
  return eslOK;

 ERROR:
  *ret_n = 0;
  return status;

}

/* Function:  cm_dsq_MPISend()
 *
 * Incept:    EPN, Tue Aug 28 15:20:16 2007
 *            
 * Purpose:   Sends a digitized sequence and it's length
 *            as a work unit to MPI process <dest> (<dest> ranges from <0..nproc-1>),
 *            tagging the message with MPI tag <tag> for MPI communicator
 *            <comm>. The receiver uses <cm_dsq_MPIRecv()> to receive the dsq.
 *            
 *            Work units are prefixed by a status code. If <dsq> is
 *            <non-NULL>, the work unit is an <eslOK> code followed by
 *            the packed MSA. If <dsq> is NULL, the work unit is an
 *            <eslEOD> code, which <cm_dsq_MPIRecv()> knows how
 *            to interpret; this is typically used for an end-of-data
 *            signal to cleanly shut down worker processes.
 *
 *            In order to minimize alloc/free cycles, caller passes a
 *            pointer to a working buffer <*buf> of size <*nalloc>
 *            characters. If necessary (i.e. if <dsq> is too big to
 *            fit), <*dsq> will be reallocated and <*nalloc> increased
 *            to the new size. As a special case, if <*dsq> is <NULL>
 *            and <*nalloc> is 0, the buffer will be allocated
 *            appropriately, but the caller is still responsible for
 *            free'ing it.
 *            
 * Args:      dsq    - digitized seq to send
 *            L      - length of dsq we're sending (dsq could extend further)
 *            dest   - MPI destination (0..nproc-1)
 *            tag    - MPI tag
 *            buf    - pointer to a working buffer 
 *            nalloc - current allocated size of <*buf>, in characters
 *
 * Returns:   <eslOK> on success; <*buf> may have been reallocated and
 *            <*nalloc> may have been increased.
 *
 * Throws:    <eslESYS> if an MPI call fails; <eslEMEM> if a malloc/realloc
 *            fails. In either case, <*buf> and <*nalloc> remain valid and useful
 *            memory (though the contents of <*buf> are undefined). 
 */
int
cm_dsq_MPISend(ESL_DSQ *dsq, int L, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc)
{
  int   status;
  int   code;
  int   sz, n, position;

  /* First, figure out the size of the work unit */
  if (MPI_Pack_size(2, MPI_INT, comm, &n) != 0) ESL_EXCEPTION(eslESYS, "mpi pack size failed"); 
  if (dsq != NULL) { 
    if (MPI_Pack_size(L+2, MPI_BYTE, comm, &sz) != 0) ESL_EXCEPTION(eslESYS, "mpi pack size failed"); 
    n += sz;
  }
  ESL_DPRINTF2(("cm_dsq_MPISend(): dsq has size %d\n", n));

  /* Make sure the buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) {
    void *tmp;
    ESL_RALLOC(*buf, tmp, sizeof(char) * n);
    *nalloc = n; 
  }
  ESL_DPRINTF2(("cm_dsq_MPISend(): buffer is ready\n"));

  /* Pack the status code, L and dsq into the buffer */
  position = 0;
  code     = (dsq == NULL) ? eslEOD : eslOK;
  if (MPI_Pack(&code, 1, MPI_INT, *buf, n, &position, comm) != 0) ESL_EXCEPTION(eslESYS, "mpi pack failed"); 
  if (dsq != NULL) {
    if (MPI_Pack(&L,  1,   MPI_INT,  *buf, n, &position, comm) != 0) ESL_EXCEPTION(eslESYS, "mpi pack failed"); 
    if (MPI_Pack(dsq, L+2, MPI_BYTE, *buf, n, &position, comm) != 0) ESL_EXCEPTION(eslESYS, "mpi pack failed"); 
  }
  ESL_DPRINTF2(("cm_dsq_MPISend(): dsq is packed into %d bytes\n", position));

  /* Send the packed profile to destination  */
  if (MPI_Send(*buf, n, MPI_PACKED, dest, tag, comm) != 0) ESL_EXCEPTION(eslESYS, "mpi send failed");
  ESL_DPRINTF2(("cm_dsq_MPISend(): dsq is sent.\n"));
  return eslOK;

 ERROR:
  return status;
}

/* Function:  cm_dsq_MPIRecv()
 *
 * Incept:    EPN, Tue Aug 28 15:29:34 2007
 *
 * Purpose:   Receives a work unit that consists of a digitized sequence
 *            and it's lenght from <source> (<0..nproc-1>, or
 *            <MPI_ANY_SOURCE>) tagged as <tag> from communicator <comm>.
 *            
 *            Work units are prefixed by a status code. If the unit's
 *            code is <eslOK> and no errors are encountered, this
 *            routine will return <eslOK> and a non-<NULL> <*ret_dsq>.
 *            If the unit's code is <eslEOD> (a shutdown signal), 
 *            this routine returns <eslEOD> and <*ret_dsq> is <NULL>.
 *            
 *            To minimize alloc/free cycles in this routine, caller
 *            passes a pointer to a buffer <*buf> of size <*nalloc>
 *            characters. These are passed by reference, because when
 *            necessary, <*buf> will be reallocated and <*nalloc>
 *            increased to the new size. As a special case, if <*buf>
 *            is <NULL> and <*nalloc> is 0, the buffer will be
 *            allocated appropriately, but the caller is still
 *            responsible for free'ing it.
 *
 *            If the packed dsq is an end-of-data signal, return
 *            <eslEOD>, and <*ret_dsq> is <NULL>.
 *            
 * Returns:   <eslOK> on success. <*ret_dsq> contains the new dsq; it
 *            is allocated here, and the caller is responsible for
 *            free'ing it.  <*buf> may have been reallocated to a
 *            larger size, and <*nalloc> may have been increased.
 *
 *
 * Throws:    <eslESYS> if an MPI call fails; <eslEMEM> if an allocation fails.
 *            In either case, <*ret_dsq> is NULL, and the <buf> and its size
 *            <*nalloc> remain valid.
 */
int
cm_dsq_MPIRecv(int source, int tag, MPI_Comm comm, char **buf, int *nalloc, ESL_DSQ **ret_dsq, int *ret_L)
{
  int         status, code;
  ESL_DSQ    *dsq  = NULL;
  int         L = 0;
  int         n;
  int         pos;
  MPI_Status  mpistatus;

  /* Probe first, because we need to know if our buffer is big enough. */
  if (MPI_Probe(source, tag, comm, &mpistatus)  != 0) ESL_XEXCEPTION(eslESYS, "mpi probe failed");
  if (MPI_Get_count(&mpistatus, MPI_PACKED, &n) != 0) ESL_XEXCEPTION(eslESYS, "mpi get count failed");

  /* Make sure the buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) {
    void *tmp;
    ESL_RALLOC(*buf, tmp, sizeof(char) * n);
    *nalloc = n; 
  }

  /* Receive the packed work unit */
  ESL_DPRINTF2(("cm_dsq_MPIRecv(): about to receive dsq.\n"));
  if (MPI_Recv(*buf, n, MPI_PACKED, source, tag, comm, &mpistatus) != 0) ESL_XEXCEPTION(eslESYS, "mpi recv failed");
  ESL_DPRINTF2(("cm_dsq_MPIRecv(): dsq has been received.\n"));

  /* Unpack it - where the first integer is a status code, OK or EOD */
  ESL_DPRINTF2(("cm_dsq_MPIRecv(): about to unpack dsq.\n"));
  pos = 0;
  if (MPI_Unpack       (*buf, n, &pos, &code,                   1, MPI_INT,           comm) != 0)  ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (code == eslEOD) { status = eslEOD; goto ERROR; }

  if (MPI_Unpack       (*buf, n, &pos, &L,                      1, MPI_INT,           comm) != 0)  ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  ESL_ALLOC(dsq, sizeof(ESL_DSQ) * (L+2));
  if (MPI_Unpack       (*buf, n, &pos, dsq,                  (L+2), MPI_BYTE,          comm) != 0)  ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  dsq[0] = dsq[(L+1)] = eslDSQ_SENTINEL; /* overwrite */
  ESL_DPRINTF2(("cm_dsq_MPIRecv(): dsq has been unpacked.\n"));

  *ret_L   = L;
  *ret_dsq = dsq;
  return eslOK;

 ERROR:
  if (dsq != NULL) free(dsq);
  *ret_dsq = NULL;
  *ret_L   = 0;
  return status;
}

/* Function:  cm_search_results_MPISend()
 *
 * Incept:    EPN, Tue Aug 28 15:51:30 2007
 *            
 * Purpose:   Send packed search results to MPI process <dest> 
 *            (<dest> ranges from <0..nproc-1>), tagging the message 
 *            with MPI tag <tag> for MPI communicator <comm>. 
 *            The receiver uses <cm_search_results_MPIRecv()> to 
 *            receive the results.
 *            
 *            Work units are prefixed by a status code. If <results> is
 *            <non-NULL>, the work unit is an <eslOK> code followed by
 *            the packed results. If <results> is NULL, the work unit is an
 *            <eslEOD> code, which <cm_search_results_MPIRecv()> knows how
 *            to interpret.
 *
 * Args:      results- search results to send
 *            dest   - MPI destination (0..nproc-1)
 *            tag    - MPI tag
 *            buf    - pointer to a working buffer 
 *            nalloc - current allocated size of <*buf>, in characters
 *
 * Returns:   <eslOK> on success; <*buf> may have been reallocated and
 *            <*nalloc> may have been increased.
 *
 * Throws:    <eslESYS> if an MPI call fails; <eslEMEM> if a malloc/realloc
 *            fails. In either case, <*buf> and <*nalloc> remain valid and useful
 *            memory (though the contents of <*buf> are undefined). 
 */
int
cm_search_results_MPISend(search_results_t *results, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc)
{
  int   status;
  int   code;
  int   sz, n, position;

  /* First, figure out the size of the work unit */
  if (MPI_Pack_size(2, MPI_INT, comm, &n) != 0) ESL_EXCEPTION(eslESYS, "mpi pack size failed"); 
  if (results != NULL) { 
    if ((status = cm_search_results_MPIPackSize(results, comm, &sz)) != eslOK) return status;
    n += sz;
  }
  ESL_DPRINTF2(("cm_search_results_MPISend(): results has size %d\n", n));

  /* Make sure the buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) {
    void *tmp;
    ESL_RALLOC(*buf, tmp, sizeof(char) * n);
    *nalloc = n; 
  }
  ESL_DPRINTF2(("cm_search_results_MPISend(): buffer is ready\n"));

  /* Pack the status code, and results into the buffer */
  position = 0;
  code     = (results == NULL) ? eslEOD : eslOK;
  if (MPI_Pack(&code, 1, MPI_INT, *buf, n, &position, comm) != 0) ESL_EXCEPTION(eslESYS, "mpi pack failed"); 
  if (results != NULL) {
    if ((status = cm_search_results_MPIPack(results, *buf, n, &position, comm)) != eslOK) return status;
  }
  ESL_DPRINTF2(("cm_search_results_MPISend(): results is packed into %d bytes\n", position));

  /* Send the packed profile to destination  */
  if (MPI_Send(*buf, n, MPI_PACKED, dest, tag, comm) != 0) ESL_EXCEPTION(eslESYS, "mpi send failed");
  ESL_DPRINTF2(("cm_search_results_MPISend(): results are sent.\n"));
  return eslOK;

 ERROR:
  return status;
}

/* Function:  cm_search_results_MPIPackSize()
 * Synopsis:  Calculates number of bytes needed to pack search
 *            results.
 * Incept:    EPN, Tue Aug 28 21:01:46 2007
 *
 * Purpose:   Calculate an upper bound on the number of bytes
 *            that <cm_search_results_MPIPack()> will need to pack
 *            search results in a packed MPI message in 
 *            communicator <comm>; return that number of bytes 
 *            in <*ret_n>. 
 *            
 *            Caller will generally use this result to determine how
 *            to allocate a buffer before starting to pack into it.
 *
 * Returns:   <eslOK> on success, and <*ret_n> contains the answer.
 *
 * Throws:    <eslESYS> if an MPI call fails, and <*ret_n> is set to 0. 
 *
 * Note:      The sizing calls here need to stay matched up with
 *            the calls in <cm_search_results_MPIPack()>.
 */
int
cm_search_results_MPIPackSize(const search_results_t *results, MPI_Comm comm, int *ret_n)
{
  int status;
  int sz;
  int n = 0;
  int i;

  status = MPI_Pack_size (1, MPI_INT, comm, &sz); n += sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  for(i = 0; i < results->num_results; i++) {
    if ((status = cm_search_result_node_MPIPackSize(&(results->data[i]), comm, &sz))  != eslOK) goto ERROR; n += sz;
  }

  *ret_n = n;
  return eslOK;

 ERROR:
  *ret_n = 0;
  return status;
}

/* Function:  cm_search_results_MPIPack()
 * Synopsis:  Packs search results into MPI buffer.
 * Incept:    EPN, Tue Aug 28 21:06:58 2007
 *
 * Purpose:   Packs search results in <results> into an 
 *            MPI packed message buffer <buf> of length <n> bytes, 
 *            starting at byte position
 *            <*position>, for MPI communicator <comm>.
 *
 * Returns:   <eslOK> on success; <buf> now contains the
 *            packed <results>, and <*position> is set to the byte
 *            immediately following the last byte of the results
 *            in <buf>. 
 *
 * Throws:    <eslESYS> if an MPI call fails; or <eslEMEM> if the
 *            buffer's length <n> is overflowed by trying to pack
 *            <results> into <buf>. In either case, the state of
 *            <buf> and <*position> is undefined, and both should
 *            be considered to be corrupted.
 *
 */
int
cm_search_results_MPIPack(const search_results_t *results, char *buf, int n, int *position, MPI_Comm comm)
{
  int status;
  int i;

  ESL_DPRINTF2(("cm_search_results_MPIPack(): ready.\n"));

  status = MPI_Pack((int *) &(results->num_results),   1, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  for (i = 0; i < results->num_results; i++) {
    status = cm_search_result_node_MPIPack(&(results->data[i]), buf, n, position, comm);  if (status != eslOK) return status;
  }
  ESL_DPRINTF2(("cm_search_results_MPIPack(): done. Packed %d bytes into buffer of size %d\n", *position, n));
  
  if (*position > n) ESL_EXCEPTION(eslEMEM, "buffer overflow");
  return eslOK;
}

/* Function:  cm_search_results_MPIUnpack()
 * Synopsis:  Unpacks search results from an MPI buffer.
 * Incept:    EPN, Wed Aug 29 05:10:20 2007
 *
 * Purpose:   Unpack a newly allocated search_results from MPI packed buffer
 *            <buf>, starting from position <*pos>, where the total length
 *            of the buffer in bytes is <n>. 
 *
 * Returns:   <eslOK> on success. <*pos> is updated to the position of
 *            the next element in <buf> to unpack (if any). <*ret_results>
 *            contains a newly allocated results, which the caller is 
 *            responsible for free'ing.
 *            
 * Throws:    <eslESYS> on an MPI call failure. <eslEMEM> on allocation failure.
 *            In either case, <*ret_results> is <NULL>, and the state of <buf>
 *            and <*pos> is undefined and should be considered to be corrupted.
 */
int
cm_search_results_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, search_results_t **ret_results)
{
  int         status;
  search_results_t *results = NULL;
  int         num_results;
  int         i;

  status = MPI_Unpack (buf, n, pos, &num_results, 1, MPI_INT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if(num_results == 0) { status = eslOK; goto ERROR; }

  results = CreateResults(num_results);
  ESL_DPRINTF2(("cm_search_results_MPIUnpack(): %d results.\n", num_results));
  for (i = 0; i < num_results; i++) {
    status = cm_search_result_node_MPIUnpack(buf, n, pos, comm, &(results->data[i]));  if (status != eslOK) return status;
    if(results->data[i].tr != NULL) ESL_DASSERT1((results->data[i].tr->n > 0));
    results->num_results++;
  }

  *ret_results = results;
  return eslOK;

 ERROR:
  if (results != NULL) FreeResults(results);
  *ret_results = NULL;
  return status;
}

/* Function:  cm_search_result_node_MPIPackSize()
 * Synopsis:  Calculates number of bytes needed to pack a search
 *            result node.
 * Incept:    EPN, Wed Aug 29 05:19:33 2007
 *
 * Purpose:   Calculate an upper bound on the number of bytes
 *            that <cm_search_result_node_MPIPack()> will need to pack
 *            a search result node in a packed MPI message in 
 *            communicator <comm>; return that number of bytes 
 *            in <*ret_n>. 
 *            
 *            Caller will generally use this result to determine how
 *            to allocate a buffer before starting to pack into it.
 *
 * Returns:   <eslOK> on success, and <*ret_n> contains the answer.
 *
 * Throws:    <eslESYS> if an MPI call fails, and <*ret_n> is set to 0. 
 *
 * Note:      The sizing calls here need to stay matched up with
 *            the calls in <cm_search_result_node_MPIPack()>.
 */
int
cm_search_result_node_MPIPackSize(const search_result_node_t *rnode, MPI_Comm comm, int *ret_n)
{
  int status;
  int sz;
  int n = 0;

  status = MPI_Pack_size (1, MPI_INT, comm, &sz);   n += 3*sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /* rnode->start, rnode->stop, rnode->bestr */
  status = MPI_Pack_size (1, MPI_FLOAT, comm, &sz); n +=   sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /* rnode->score */
  status = MPI_Pack_size (1, MPI_INT, comm, &sz);   n += 2*sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /* has_tr, has_pcodes */
  if(rnode->tr != NULL) {
    if ((status = cm_parsetree_MPIPackSize(rnode->tr, comm, &sz))  != eslOK) goto ERROR; n += sz;
  }
  if(rnode->pcode != NULL) { 
    if ((status = esl_mpi_PackOptSize(rnode->pcode, -1, MPI_CHAR, comm, &sz)) != eslOK) goto ERROR; n += sz;
  }

  *ret_n = n;
  return eslOK;

 ERROR:
  *ret_n = 0;
  return status;
}

/* Function:  cm_search_result_node_MPIPack()
 * Synopsis:  Packs search result node into MPI buffer.
 * Incept:    EPN, Wed Aug 29 05:22:32 2007
 *
 * Purpose:   Packs search result node in <rnode> into an 
 *            MPI packed message buffer <buf> of length <n> bytes, 
 *            starting at byte position
 *            <*position>, for MPI communicator <comm>.
 *
 * Returns:   <eslOK> on success; <buf> now contains the
 *            packed <rnode>, and <*position> is set to the byte
 *            immediately following the last byte of the results
 *            in <buf>. 
 *
 * Throws:    <eslESYS> if an MPI call fails; or <eslEMEM> if the
 *            buffer's length <n> is overflowed by trying to pack
 *            <rnode> into <buf>. In either case, the state of
 *            <buf> and <*position> is undefined, and both should
 *            be considered to be corrupted.
 *
 */
int
cm_search_result_node_MPIPack(const search_result_node_t *rnode, char *buf, int n, int *position, MPI_Comm comm)
{
  int status;
  int has_tr;
  int has_pcodes;

  ESL_DPRINTF2(("cm_search_result_node_MPIPack(): ready.\n"));
  
  status = MPI_Pack((int *) &(rnode->start),   1, MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *) &(rnode->stop),    1, MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *) &(rnode->bestr),   1, MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *) &(rnode->score),   1, MPI_FLOAT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");

  has_tr     = (rnode->tr != NULL) ? TRUE : FALSE;
  has_pcodes = (rnode->pcode != NULL) ? TRUE : FALSE;
  status = MPI_Pack(&has_tr,                 1, MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(&has_pcodes,             1, MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");

  if(has_tr) { 
    status = cm_parsetree_MPIPack((Parsetree_t *) rnode->tr, buf, n, position, comm);  if (status != eslOK) return status;
  }
  if(has_pcodes) {
    /* we call PackOpt, even though we know we should have valid posterior codes */
    status = esl_mpi_PackOpt(rnode->pcode, -1, MPI_CHAR, buf, n, position,  comm); if (status != eslOK) return status;
  }

  ESL_DPRINTF2(("cm_search_result_node_MPIPack(): done. Packed %d bytes into buffer of size %d\n", *position, n));

  if (*position > n) ESL_EXCEPTION(eslEMEM, "buffer overflow");
  return eslOK;
}

/* Function:  cm_search_result_node_MPIUnpack()
 * Synopsis:  Unpacks search result node from an MPI buffer.
 * Incept:    EPN, Wed Aug 29 05:10:20 2007
 *
 * Purpose:   Unpack a newly allocated search result node from MPI packed buffer
 *            <buf>, starting from position <*pos>, where the total length
 *            of the buffer in bytes is <n>. 
 *
 * Returns:   <eslOK> on success. <*pos> is updated to the position of
 *            the next element in <buf> to unpack (if any). <*ret_rnode>
 *            contains a newly allocated result node, which the caller is 
 *            responsible for free'ing.
 *            
 * Throws:    <eslESYS> on an MPI call failure. <eslEMEM> on allocation failure.
 *            In either case, <*ret_rnode> is <NULL>, and the state of <buf>
 *            and <*pos> is undefined and should be considered to be corrupted.
 */
int
cm_search_result_node_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, search_result_node_t *ret_rnode)
{
  int status;
  search_result_node_t rnode;
  int start, stop, bestr;
  float score;
  int has_tr = FALSE;
  int has_pcodes = FALSE;
  Parsetree_t *tr = NULL;
  char *pcode = NULL;

  status = MPI_Unpack (buf, n, pos, &start, 1, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &stop,  1, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &bestr, 1, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &score, 1, MPI_FLOAT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &has_tr, 1, MPI_INT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &has_pcodes, 1, MPI_INT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  rnode.start = start;
  rnode.stop  = stop;
  rnode.bestr = bestr;
  rnode.score = score;
  rnode.tr    = NULL;
  rnode.pcode = NULL;

  /* optionally, unpack a parsetree */
  if(has_tr) {
    status   = cm_parsetree_MPIUnpack(buf, n, pos, comm, &tr);  if (status != eslOK) return status;
    rnode.tr = tr;
  }
  if(has_pcodes) {
    status = esl_mpi_UnpackOpt(buf, n, pos, (void **) &(pcode), NULL, MPI_CHAR, comm); if (status != eslOK) goto ERROR;;
    rnode.pcode = pcode;
  }
  
  *ret_rnode = rnode;
  return eslOK;
  
 ERROR:
  if(tr != NULL) FreeParsetree(tr);
  if(pcode != NULL) free(pcode);
  ret_rnode = NULL;
  return status;
}

/* Function:  cm_seqs_to_aln_MPISend()
 *
 * Incept:    EPN, Mon Sep  3 14:58:12 2007
 *            
 * Purpose:   Send packed seqs_to_aln to MPI process <dest> 
 *            (<dest> ranges from <0..nproc-1>), tagging the message 
 *            with MPI tag <tag> for MPI communicator <comm>. 
 *            The receiver uses <cm_seqs_to_aln_MPIRecv()> to 
 *            receive the results.
 *            
 *            Work units are prefixed by a status code. If <results> is
 *            <non-NULL>, the work unit is an <eslOK> code followed by
 *            the packed results. If <results> is NULL, the work unit is an
 *            <eslEOD> code, which <cm_seqs_to_aln_MPIRecv()> knows how
 *            to interpret.
 *
 * Args:      seqs_to_aln  - seqs_to_aln_t object to send
 *            offset       - index of first to send (send seqs_to_aln[offset] first)
 *            nseq_to_send - number to send
 *            dest         - MPI destination (0..nproc-1)
 *            tag          - MPI tag
 *            buf          - pointer to a working buffer 
 *            nalloc       - current allocated size of <*buf>, in characters
 *
 * Returns:   <eslOK> on success; <*buf> may have been reallocated and
 *            <*nalloc> may have been increased.
 *
 * Throws:    <eslESYS> if an MPI call fails; <eslEMEM> if a malloc/realloc
 *            fails. In either case, <*buf> and <*nalloc> remain valid and useful
 *            memory (though the contents of <*buf> are undefined). 
 */
int
cm_seqs_to_aln_MPISend(seqs_to_aln_t *seqs_to_aln, int offset, int nseq_to_send, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc)
{
  int   status;
  int   code;
  int   sz, n, position;

  /* First, figure out the size of the work unit */
  if (MPI_Pack_size(2, MPI_INT, comm, &n) != 0) ESL_EXCEPTION(eslESYS, "mpi pack size failed"); 
  /* this is size of the work unit, and the status code */

  if (seqs_to_aln != NULL) { 
    if ((status = cm_seqs_to_aln_MPIPackSize(seqs_to_aln, offset, nseq_to_send, comm, &sz)) != eslOK) return status;
    n += sz;
  }
  ESL_DPRINTF2(("cm_seqs_to_aln_MPISend(): seqs_to_aln has size %d\n", n));

  /* Make sure the buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) {
    void *tmp;
    ESL_RALLOC(*buf, tmp, sizeof(char) * n);
    *nalloc = n; 
  }
  ESL_DPRINTF2(("cm_seqs_to_aln_MPISend(): buffer is ready\n"));

  /* Pack the status code, and seqs_to_aln into the buffer */
  position = 0;
  code     = (seqs_to_aln == NULL) ? eslEOD : eslOK;
  if (MPI_Pack(&code, 1, MPI_INT, *buf, n, &position, comm) != 0) ESL_EXCEPTION(eslESYS, "mpi pack failed"); 
  if (seqs_to_aln != NULL) {
    if ((status = cm_seqs_to_aln_MPIPack(seqs_to_aln, offset, nseq_to_send, *buf, n, &position, comm)) != eslOK) return status;
  }
  ESL_DPRINTF2(("cm_seqs_to_aln_MPISend(): seqs_to_aln is packed into %d bytes\n", position));

  /* Send the packed profile to destination  */
  if (MPI_Send(*buf, n, MPI_PACKED, dest, tag, comm) != 0) ESL_EXCEPTION(eslESYS, "mpi send failed");
  ESL_DPRINTF2(("cm_seqs_to_aln_MPISend(): seqs_to_aln are sent.\n"));
  return eslOK;

 ERROR:
  return status;
}

/* Function:  cm_seqs_to_aln_MPIRecv()
 *
 * Incept:    EPN, Mon Sep  3 19:42:47 2007
 *
 * Purpose:   Receives a work unit that consists of a seqs_to_aln data structure.
 *            from <source> (<0..nproc-1>, or <MPI_ANY_SOURCE>) tagged as <tag> 
 *            from communicator <comm>.
 *            
 *            Work units are prefixed by a status code. If the unit's
 *            code is <eslOK> and no errors are encountered, this
 *            routine will return <eslOK> and a non-<NULL> <*ret_seqs_to_aln>.
 *            If the unit's code is <eslEOD> (a shutdown signal), 
 *            this routine returns <eslEOD> and <*ret_seqs_to_aln> is <NULL>.
 *            
 *            To minimize alloc/free cycles in this routine, caller
 *            passes a pointer to a buffer <*buf> of size <*nalloc>
 *            characters. These are passed by reference, because when
 *            necessary, <*buf> will be reallocated and <*nalloc>
 *            increased to the new size. As a special case, if <*buf>
 *            is <NULL> and <*nalloc> is 0, the buffer will be
 *            allocated appropriately, but the caller is still
 *            responsible for free'ing it.
 *
 *            If the packed seqs_to_aln is an end-of-data signal, return
 *            <eslEOD>, and <*ret_seqs_to_aln> is <NULL>.
 *            
 * Returns:   <eslOK> on success. <*ret_seqs_to_aln> contains the new 
 *            seqs_to_aln object; it is allocated here, and the caller 
 *            is responsible for free'ing it.  <*buf> may have been 
 *            reallocated to a larger size, and <*nalloc> may have been 
 *            increased.
 *
 * Throws:    <eslESYS> if an MPI call fails; <eslEMEM> if an allocation fails.
 *            In either case, <*ret_seqs_to_aln> is NULL, and the <buf> and its size
 *            <*nalloc> remain valid.
 */
int
cm_seqs_to_aln_MPIRecv(const ESL_ALPHABET *abc, int source, int tag, MPI_Comm comm, char **buf, int *nalloc, seqs_to_aln_t **ret_seqs_to_aln)
{
  int         status, code;
  seqs_to_aln_t *seqs_to_aln = NULL;
  int         n;
  int         pos;
  MPI_Status  mpistatus;

  /* Probe first, because we need to know if our buffer is big enough. */
  if (MPI_Probe(source, tag, comm, &mpistatus)  != 0) ESL_XEXCEPTION(eslESYS, "mpi probe failed");
  if (MPI_Get_count(&mpistatus, MPI_PACKED, &n) != 0) ESL_XEXCEPTION(eslESYS, "mpi get count failed");

  /* Make sure the buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) {
    void *tmp;
    ESL_RALLOC(*buf, tmp, sizeof(char) * n);
    *nalloc = n; 
  }

  /* Receive the packed work unit */
  ESL_DPRINTF2(("cm_seqs_to_aln_MPIRecv(): about to receive seqs_to_aln.\n"));
  if (MPI_Recv(*buf, n, MPI_PACKED, source, tag, comm, &mpistatus) != 0) ESL_XEXCEPTION(eslESYS, "mpi recv failed");
  ESL_DPRINTF2(("cm_seqs_to_aln_MPIRecv(): seqs_to_aln has been received.\n"));

  /* Unpack it - where the first integer is a status code, OK or EOD */
  ESL_DPRINTF2(("cm_seqs_to_aln_MPIRecv(): about to unpack seqs_to_aln.\n"));
  pos = 0;
  if (MPI_Unpack       (*buf, n, &pos, &code,                   1, MPI_INT,           comm) != 0)  ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (code == eslEOD) { status = eslEOD; goto ERROR; }

  return cm_seqs_to_aln_MPIUnpack(abc, *buf, *nalloc, &pos, comm, ret_seqs_to_aln);

 ERROR:
  if (seqs_to_aln != NULL) FreeSeqsToAln(seqs_to_aln);
  *ret_seqs_to_aln = NULL;
  return status;
}

/* Function:  cm_seqs_to_aln_MPIPackSize()
 * Synopsis:  Calculates number of bytes needed to pack a
 *            seqs_to_aln object.
 * Incept:    EPN, Mon Sep  3 19:49:31 2007
 *
 * Purpose:   Calculate an upper bound on the number of bytes
 *            that <cm_seqs_to_aln_MPIPack()> will need to pack
 *            search results in a packed MPI message in 
 *            communicator <comm>; return that number of bytes 
 *            in <*ret_n>. 
 *            
 *            Caller will generally use this result to determine how
 *            to allocate a buffer before starting to pack into it.
 *
 * Returns:   <eslOK> on success, and <*ret_n> contains the answer.
 *
 * Throws:    <eslESYS> if an MPI call fails, and <*ret_n> is set to 0. 
 *
 * Note:      The sizing calls here need to stay matched up with
 *            the calls in <cm_seqs_to_aln_MPIPack()>.
 */
int
cm_seqs_to_aln_MPIPackSize(const seqs_to_aln_t *seqs_to_aln, int offset, int nseq_to_pack, MPI_Comm comm, int *ret_n)
{
  int status;
  int sz;
  int n = 0;
  int i;

  /* determine what information we have in the seqs_to_aln object */
  int has_sq      = TRUE;
  int has_tr      = TRUE;
  int has_cp9_tr  = TRUE;
  int has_post    = TRUE;  
  int has_sc      = TRUE;  
  int has_pp      = TRUE;  
  int has_struct_sc = TRUE;  

  /* careful, individual sq, tr, cp9_tr, postcode ptrs may be NULL even if ptr to ptrs is non-NULL,
   * example sq != NULL and sq[i] == NULL is possible. */

  if(seqs_to_aln->sq == NULL) has_sq = FALSE;
  else
    for (i = offset; i < offset + nseq_to_pack; i++) 
      if(seqs_to_aln->sq[i] == NULL) { has_sq = FALSE; break; }

  if(seqs_to_aln->tr == NULL) has_tr = FALSE;
  else
    for (i = offset; i < offset + nseq_to_pack; i++) 
      if(seqs_to_aln->tr[i] == NULL) { has_tr = FALSE; break; }

  if(seqs_to_aln->cp9_tr == NULL) has_cp9_tr = FALSE;
  else
    for (i = offset; i < offset + nseq_to_pack; i++) 
      if(seqs_to_aln->cp9_tr[i] == NULL) { has_cp9_tr = FALSE; break; }

  if(seqs_to_aln->postcode == NULL) has_post = FALSE;
  else
    for (i = offset; i < offset + nseq_to_pack; i++) 
      if(seqs_to_aln->postcode[i] == NULL) { has_post = FALSE; break; }

  if(seqs_to_aln->sc == NULL) has_sc = FALSE;
  else
    for (i = offset; i < offset + nseq_to_pack; i++) 
      if(! NOT_IMPOSSIBLE(seqs_to_aln->sc[i])) { has_sc = FALSE; break; }

  if(seqs_to_aln->pp == NULL) has_pp = FALSE;
  else
    for (i = offset; i < offset + nseq_to_pack; i++) 
      if(! NOT_IMPOSSIBLE(seqs_to_aln->pp[i])) { has_pp = FALSE; break; }

  if(seqs_to_aln->struct_sc == NULL) has_struct_sc = FALSE;
  else
    for (i = offset; i < offset + nseq_to_pack; i++) { 
      if(! NOT_IMPOSSIBLE(seqs_to_aln->struct_sc[i])) { has_struct_sc = FALSE; break; }
    }
  status = MPI_Pack_size (1, MPI_INT, comm, &sz); n += 8*sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /* nseq_to_pack, has_sq, has_tr, has_cp9_tr, has_post, has_sc, has_pp, has_struct_sc */

  if(has_sq) {
    for(i = offset; i < offset + nseq_to_pack; i++) {
      if ((status = cm_digitized_sq_MPIPackSize(seqs_to_aln->sq[i], comm, &sz))  != eslOK) goto ERROR; n += sz;
    }
  }
  if(has_tr) {
    for(i = offset; i < offset + nseq_to_pack; i++) {
      if ((status = cm_parsetree_MPIPackSize(seqs_to_aln->tr[i], comm, &sz))  != eslOK) goto ERROR; n += sz;
    }
  }
  if(has_cp9_tr) {
    for(i = offset; i < offset + nseq_to_pack; i++) {
      if ((status = cm_cp9trace_MPIPackSize(seqs_to_aln->cp9_tr[i], comm, &sz))  != eslOK) goto ERROR; n += sz;
    }
  }
  if(has_post) {
    for(i = offset; i < offset + nseq_to_pack; i++) {
      if ((status = esl_mpi_PackOptSize(seqs_to_aln->postcode[i], -1, MPI_CHAR, comm, &sz)) != eslOK) goto ERROR; n += sz;
    }
  }
  if(has_sc) 
    if ((status = MPI_Pack_size(nseq_to_pack, MPI_FLOAT, comm, &sz)) != eslOK) goto ERROR; n += sz;
  if(has_pp) 
    if ((status = MPI_Pack_size(nseq_to_pack, MPI_FLOAT, comm, &sz)) != eslOK) goto ERROR; n += sz;
  if(has_struct_sc) 
    if ((status = MPI_Pack_size(nseq_to_pack, MPI_FLOAT, comm, &sz)) != eslOK) goto ERROR; n += sz;

  *ret_n = n;
  return eslOK;

 ERROR:
  *ret_n = 0;
  return status;
}

/* Function:  cm_seqs_to_aln_MPIPack()
 * Synopsis:  Packs seqs_to_aln into MPI buffer.
 * Incept:    EPN, Mon Sep  3 15:04:33 2007
 *
 * Purpose:   Packs some or all of the sequences in <seqs_to_aln> into an 
 *            MPI packed message buffer <buf> of length <n> bytes, 
 *            starting at byte position
 *            <*position>, for MPI communicator <comm>.
 *
 * Returns:   <eslOK> on success; <buf> now contains the
 *            packed <results>, and <*position> is set to the byte
 *            immediately following the last byte of the results
 *            in <buf>. 
 *
 * Throws:    <eslESYS> if an MPI call fails; or <eslEMEM> if the
 *            buffer's length <n> is overflowed by trying to pack
 *            <results> into <buf>. In either case, the state of
 *            <buf> and <*position> is undefined, and both should
 *            be considered to be corrupted.
 *
 */
int
cm_seqs_to_aln_MPIPack(const seqs_to_aln_t *seqs_to_aln, int offset, int nseq_to_pack, char *buf, int n, int *position, MPI_Comm comm)
{
  int status;
  int i;

  ESL_DPRINTF2(("cm_seqs_to_aln_MPIPack(): ready.\n"));
  ESL_DASSERT1((seqs_to_aln->nseq >= (offset + nseq_to_pack)));

  /* determine what information we have in the seqs_to_aln object */
  int has_sq      = TRUE;
  int has_tr      = TRUE;
  int has_cp9_tr  = TRUE;
  int has_post    = TRUE;  
  int has_sc      = TRUE;  
  int has_pp      = TRUE;  
  int has_struct_sc = TRUE;  
  /* careful, individual sq, tr, cp9_tr, postcode ptrs may be NULL even if ptr to ptrs is non-NULL,
   * example sq != NULL and sq[i] == NULL is possible. */

  if(seqs_to_aln->sq == NULL) has_sq = FALSE;
  else
    for (i = offset; i < offset + nseq_to_pack; i++) 
      if(seqs_to_aln->sq[i] == NULL) { has_sq = FALSE; break; }

  if(seqs_to_aln->tr == NULL) has_tr = FALSE;
  else
    for (i = offset; i < offset + nseq_to_pack; i++) 
      if(seqs_to_aln->tr[i] == NULL) { has_tr = FALSE; break; }

  if(seqs_to_aln->cp9_tr == NULL) has_cp9_tr = FALSE;
  else
    for (i = offset; i < offset + nseq_to_pack; i++) 
      if(seqs_to_aln->cp9_tr[i] == NULL) { has_cp9_tr = FALSE; break; }

  if(seqs_to_aln->postcode == NULL) has_post = FALSE;
  else
    for (i = offset; i < offset + nseq_to_pack; i++) 
      if(seqs_to_aln->postcode[i] == NULL) { has_post = FALSE; break; }

  if(seqs_to_aln->sc == NULL) has_sc = FALSE;
  else
    for (i = offset; i < offset + nseq_to_pack; i++) 
      if(! NOT_IMPOSSIBLE(seqs_to_aln->sc[i])) { has_sc = FALSE; break; }

  if(seqs_to_aln->pp == NULL) has_pp = FALSE;
  else
    for (i = offset; i < offset + nseq_to_pack; i++) 
      if(! NOT_IMPOSSIBLE(seqs_to_aln->pp[i])) { has_pp = FALSE; break; }

  if(seqs_to_aln->struct_sc == NULL) has_struct_sc = FALSE;
  else
    for (i = offset; i < offset + nseq_to_pack; i++) 
      if(! NOT_IMPOSSIBLE(seqs_to_aln->struct_sc[i])) { has_struct_sc = FALSE; break; }


  status = MPI_Pack((int *) &(nseq_to_pack), 1, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *) &(has_sq),       1, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *) &(has_tr),       1, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *) &(has_cp9_tr),   1, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *) &(has_post),     1, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *) &(has_sc),       1, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *) &(has_pp),       1, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *) &(has_struct_sc),1, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");

  if(has_sq)
    for (i = offset; i < offset + nseq_to_pack; i++) {
      status = cm_digitized_sq_MPIPack(seqs_to_aln->sq[i], buf, n, position, comm);  if (status != eslOK) return status;
    }

  if(has_tr)
    for (i = offset; i < offset + nseq_to_pack; i++) {
      status = cm_parsetree_MPIPack(seqs_to_aln->tr[i], buf, n, position, comm);  if (status != eslOK) return status;
    }

  if(has_cp9_tr)
    for (i = offset; i < offset + nseq_to_pack; i++) {
      status = cm_cp9trace_MPIPack(seqs_to_aln->cp9_tr[i], buf, n, position, comm);  if (status != eslOK) return status;
    }

  if(has_post)
    for (i = offset; i < offset + nseq_to_pack; i++) {
      /* we call PackOpt, even though we know we should have valid posterior codes */
      status = esl_mpi_PackOpt(seqs_to_aln->postcode[i], -1, MPI_CHAR, buf, n, position,  comm); if (status != eslOK) return status;
    }

  if(has_sc)
    status = MPI_Pack((seqs_to_aln->sc + (offset * sizeof(float))), nseq_to_pack, MPI_FLOAT, buf, n, position, comm); if (status != eslOK) return status;

  if(has_pp)
    status = MPI_Pack((seqs_to_aln->pp + (offset * sizeof(float))), nseq_to_pack, MPI_FLOAT, buf, n, position, comm); if (status != eslOK) return status;

  if(has_struct_sc)
    status = MPI_Pack((seqs_to_aln->struct_sc + (offset * sizeof(float))), nseq_to_pack, MPI_FLOAT, buf, n, position, comm); if (status != eslOK) return status;

  ESL_DPRINTF2(("cm_seqs_to_aln_MPIPack(): done. Packed %d bytes into buffer of size %d\n", *position, n));
  
  if (*position > n) ESL_EXCEPTION(eslEMEM, "buffer overflow");
  return eslOK;
}

/* Function:  cm_seqs_to_aln_MPIUnpack()
 * Synopsis:  Unpacks sequences seqs_to_aln from an MPI buffer.
 * Incept:    EPN, Mon Sep  3 15:48:11 2007
 *
 * Purpose:   Unpack a newly allocated seqs_to_aln from MPI packed buffer
 *            <buf>, starting from position <*pos>, where the total length
 *            of the buffer in bytes is <n>. 
 *
 * Returns:   <eslOK> on success. <*pos> is updated to the position of
 *            the next element in <buf> to unpack (if any). <*ret_seqs_to_aln>
 *            contains a newly allocated seqs_to_aln, which the caller is 
 *            responsible for free'ing.
 *            
 * Throws:    <eslESYS> on an MPI call failure. <eslEMEM> on allocation failure.
 *            In either case, <*ret_seqs_to_aln> is <NULL>, and the state of <buf>
 *            and <*pos> is undefined and should be considered to be corrupted.
 */
int
cm_seqs_to_aln_MPIUnpack(const ESL_ALPHABET *abc, char *buf, int n, int *pos, MPI_Comm comm, seqs_to_aln_t **ret_seqs_to_aln)
{
  int         status;
  seqs_to_aln_t *seqs_to_aln = NULL;
  int         num_seqs_to_aln;
  int         i;
  int         has_sq;
  int         has_tr;
  int         has_cp9_tr;
  int         has_post;
  int         has_sc;
  int         has_pp;
  int         has_struct_sc;

  status = MPI_Unpack (buf, n, pos, &num_seqs_to_aln, 1, MPI_INT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &has_sq,          1, MPI_INT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &has_tr,          1, MPI_INT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &has_cp9_tr,      1, MPI_INT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &has_post,        1, MPI_INT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &has_sc,          1, MPI_INT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &has_pp,          1, MPI_INT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &has_struct_sc,   1, MPI_INT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if(num_seqs_to_aln == 0) { status = eslOK; goto ERROR; }

  seqs_to_aln = CreateSeqsToAln(num_seqs_to_aln, FALSE);
  ESL_DPRINTF2(("cm_seqs_to_aln_MPIUnpack(): %d seqs_to_aln.\n", num_seqs_to_aln));

  if(has_sq)  
    for (i = 0; i < num_seqs_to_aln; i++) {
      status = cm_digitized_sq_MPIUnpack(abc, buf, n, pos, comm, &(seqs_to_aln->sq[i]));  if (status != eslOK) goto ERROR;;
    }
  else if(seqs_to_aln->sq != NULL) { free(seqs_to_aln->sq); seqs_to_aln->sq = NULL; }

  if(has_tr) {
    ESL_ALLOC(seqs_to_aln->tr, sizeof(Parsetree_t *) * num_seqs_to_aln);
    for (i = 0; i < num_seqs_to_aln; i++) {
      status = cm_parsetree_MPIUnpack(buf, n, pos, comm, &(seqs_to_aln->tr[i]));  if (status != eslOK) goto ERROR;;
    }
  }

  if(has_cp9_tr) {
    ESL_ALLOC(seqs_to_aln->cp9_tr, sizeof(CP9trace_t *) * num_seqs_to_aln);
    for (i = 0; i < num_seqs_to_aln; i++) {
      status = cm_cp9trace_MPIUnpack(buf, n, pos, comm, &(seqs_to_aln->cp9_tr[i]));  if (status != eslOK) goto ERROR;;
    }
  }

  if(has_post) {
    ESL_ALLOC(seqs_to_aln->postcode, sizeof(char *) * num_seqs_to_aln);
    for (i = 0; i < num_seqs_to_aln; i++) {
      status = esl_mpi_UnpackOpt(buf, n, pos, (void **) &(seqs_to_aln->postcode[i]), NULL, MPI_CHAR, comm); if (status != eslOK) goto ERROR;;
    }
  }

  if(has_sc) {
    ESL_ALLOC(seqs_to_aln->sc, sizeof(float) * num_seqs_to_aln);
    status = MPI_Unpack(buf, n, pos, seqs_to_aln->sc, num_seqs_to_aln, MPI_FLOAT, comm); if (status != eslOK) goto ERROR;;
  }

  if(has_pp) {
    ESL_ALLOC(seqs_to_aln->pp, sizeof(float) * num_seqs_to_aln);
    status = MPI_Unpack(buf, n, pos, seqs_to_aln->pp, num_seqs_to_aln, MPI_FLOAT, comm); if (status != eslOK) goto ERROR;;
  }

  if(has_struct_sc) {
    ESL_ALLOC(seqs_to_aln->struct_sc, sizeof(float) * num_seqs_to_aln);
    status = MPI_Unpack(buf, n, pos, seqs_to_aln->struct_sc, num_seqs_to_aln, MPI_FLOAT, comm); if (status != eslOK) goto ERROR;;
  }

  seqs_to_aln->nseq = num_seqs_to_aln;
  *ret_seqs_to_aln = seqs_to_aln;
  return eslOK;

 ERROR:
  if (seqs_to_aln != NULL) FreeSeqsToAln(seqs_to_aln);
  *ret_seqs_to_aln = NULL;
  return status;
}

/* Function:  cm_parsetree_MPIPackSize()
 * Synopsis:  Calculates number of bytes needed to pack a 
 *            parsetree.
 * Incept:    EPN, Wed Aug 29 05:44:28 2007
 *
 * Purpose:   Calculate an upper bound on the number of bytes
 *            that <cm_parsetree_MPIPack()> will need to pack
 *            a parsetree in a packed MPI message in 
 *            communicator <comm>; return that number of bytes 
 *            in <*ret_n>. 
 *            
 *            Caller will generally use this result to determine how
 *            to allocate a buffer before starting to pack into it.
 *
 * Returns:   <eslOK> on success, and <*ret_n> contains the answer.
 *
 * Throws:    <eslESYS> if an MPI call fails, and <*ret_n> is set to 0. 
 *
 * Note:      The sizing calls here need to stay matched up with
 *            the calls in <cm_parsetree_MPIPack()>.
 */
int
cm_parsetree_MPIPackSize(const Parsetree_t *tr, MPI_Comm comm, int *ret_n)
{
  int status;
  int sz;
  int n = 0;

  status = MPI_Pack_size(1,     MPI_INT, comm, &sz); n += 2*sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  status = MPI_Pack_size(tr->n, MPI_INT, comm, &sz); n += 7*sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");

  *ret_n = n;
  return eslOK;

 ERROR:
  *ret_n = 0;
  return status;
}

/* Function:  cm_parsetree_MPIPack()
 * Synopsis:  Packs parsetree into MPI buffer.
 * Incept:    EPN, Wed Aug 29 05:22:32 2007
 *
 * Purpose:   Packs parsetree <tr> into an 
 *            MPI packed message buffer <buf> of length <n> bytes, 
 *            starting at byte position
 *            <*position>, for MPI communicator <comm>.
 *
 * Returns:   <eslOK> on success; <buf> now contains the
 *            packed <tr>, and <*position> is set to the byte
 *            immediately following the last byte of the results
 *            in <buf>. 
 *
 * Throws:    <eslESYS> if an MPI call fails; or <eslEMEM> if the
 *            buffer's length <n> is overflowed by trying to pack
 *            <rnode> into <buf>. In either case, the state of
 *            <buf> and <*position> is undefined, and both should
 *            be considered to be corrupted.
 *
 */
int
cm_parsetree_MPIPack(const Parsetree_t *tr, char *buf, int n, int *position, MPI_Comm comm)
{
  int status;

  ESL_DPRINTF2(("cm_parsetree_MPIPack(): ready.\n"));
  
  status = MPI_Pack((int *) &(tr->n),        1, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *) &(tr->memblock), 1, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(tr->emitl,           tr->n, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(tr->emitr,           tr->n, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(tr->state,           tr->n, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(tr->mode,            tr->n, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(tr->nxtl,            tr->n, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(tr->nxtr,            tr->n, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(tr->prv,             tr->n, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");

  ESL_DPRINTF2(("cm_parsetree_MPIPack(): done. Packed %d bytes into buffer of size %d\n", *position, n));

  if (*position > n) ESL_EXCEPTION(eslEMEM, "buffer overflow");
  return eslOK;
}

/* Function:  cm_parsetree_MPIUnpack()
 * Synopsis:  Unpacks parsetree from an MPI buffer.
 * Incept:    EPN, Wed Aug 29 05:10:20 2007
 *
 * Purpose:   Unpack a newly allocated parsetree node from MPI packed buffer
 *            <buf>, starting from position <*pos>, where the total length
 *            of the buffer in bytes is <n>. 
 *
 * Returns:   <eslOK> on success. <*pos> is updated to the position of
 *            the next element in <buf> to unpack (if any). <*ret_tr>
 *            contains a newly allocated parsetree, which the caller is 
 *            responsible for free'ing.
 *            
 * Throws:    <eslESYS> on an MPI call failure. <eslEMEM> on allocation failure.
 *            In either case, <*ret_tr> is <NULL>, and the state of <buf>
 *            and <*pos> is undefined and should be considered to be corrupted.
 */
int
cm_parsetree_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, Parsetree_t **ret_tr)
{
  int status;
  Parsetree_t *tr = NULL;
  int tr_n, memblock; /* memblock is likely irrelevant, the tree probably won't grow */

  status = MPI_Unpack (buf, n, pos, &tr_n,        1, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  tr = CreateParsetree(tr_n);
  status = MPI_Unpack (buf, n, pos, &memblock,    1, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  tr->memblock = memblock;
  tr->n = tr_n;

  status = MPI_Unpack (buf, n, pos, tr->emitl, tr_n, MPI_INT,  comm); if (status != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, tr->emitr, tr_n, MPI_INT,  comm); if (status != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, tr->state, tr_n, MPI_INT,  comm); if (status != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, tr->mode,  tr_n, MPI_INT,  comm); if (status != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, tr->nxtl,  tr_n, MPI_INT,  comm); if (status != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, tr->nxtr,  tr_n, MPI_INT,  comm); if (status != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, tr->prv,   tr_n, MPI_INT,  comm); if (status != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  *ret_tr = tr;
  return eslOK;

 ERROR:
  if (tr != NULL) 
    FreeParsetree(tr);
  *ret_tr = NULL;
  return status;
}

/* Function:  cm_cp9trace_MPIPackSize()
 * Synopsis:  Calculates number of bytes needed to pack a 
 *            cp9trace.
 * Incept:    EPN, Mon Sep  3 14:49:41 2007
 *
 * Purpose:   Calculate an upper bound on the number of bytes
 *            that <cm_cp9trace_MPIPack()> will need to pack
 *            a parsetree in a packed MPI message in 
 *            communicator <comm>; return that number of bytes 
 *            in <*ret_n>. 
 *            
 *            Caller will generally use this result to determine how
 *            to allocate a buffer before starting to pack into it.
 *
 * Returns:   <eslOK> on success, and <*ret_n> contains the answer.
 *
 * Throws:    <eslESYS> if an MPI call fails, and <*ret_n> is set to 0. 
 *
 * Note:      The sizing calls here need to stay matched up with
 *            the calls in <cm_cp9trace_MPIPack()>.
 */
int
cm_cp9trace_MPIPackSize(const CP9trace_t *cp9_tr, MPI_Comm comm, int *ret_n)
{
  int status;
  int sz;
  int n = 0;

  status = MPI_Pack_size(1,            MPI_INT,  comm, &sz); n += sz;   if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  status = MPI_Pack_size(cp9_tr->tlen, MPI_INT,  comm, &sz); n += 2*sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  status = MPI_Pack_size(cp9_tr->tlen, MPI_CHAR, comm, &sz); n += sz;   if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");

  *ret_n = n;
  return eslOK;

 ERROR:
  *ret_n = 0;
  return status;
}

/* Function:  cm_cp9trace_MPIPack()
 * Synopsis:  Packs cp9trace into MPI buffer.
 * Incept:    EPN, Wed Aug 29 05:22:32 2007
 *
 * Purpose:   Packs cp9trace <cp9_tr> into an 
 *            MPI packed message buffer <buf> of length <n> bytes, 
 *            starting at byte position
 *            <*position>, for MPI communicator <comm>.
 *
 * Returns:   <eslOK> on success; <buf> now contains the
 *            packed <cp9_tr>, and <*position> is set to the byte
 *            immediately following the last byte of the results
 *            in <buf>. 
 *
 * Throws:    <eslESYS> if an MPI call fails; or <eslEMEM> if the
 *            buffer's length <n> is overflowed by trying to pack
 *            <rnode> into <buf>. In either case, the state of
 *            <buf> and <*position> is undefined, and both should
 *            be considered to be corrupted.
 *
 */
int
cm_cp9trace_MPIPack(const CP9trace_t *cp9_tr, char *buf, int n, int *position, MPI_Comm comm)
{
  int status;

  ESL_DPRINTF2(("cm_cp9trace_MPIPack(): ready.\n"));
  
  status = MPI_Pack((int *) &(cp9_tr->tlen), 1, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(cp9_tr->statetype,       cp9_tr->tlen, MPI_CHAR, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(cp9_tr->nodeidx,         cp9_tr->tlen, MPI_INT,  buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(cp9_tr->pos,             cp9_tr->tlen, MPI_INT,  buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");

  ESL_DPRINTF2(("cm_cp9trace_MPIPack(): done. Packed %d bytes into buffer of size %d\n", *position, n));

  if (*position > n) ESL_EXCEPTION(eslEMEM, "buffer overflow");
  return eslOK;
}

/* Function:  cm_cp9trace_MPIUnpack()
 * Synopsis:  Unpacks cp9trace from an MPI buffer.
 * Incept:    EPN, Wed Aug 29 05:10:20 2007
 *
 * Purpose:   Unpack a newly allocated cp9trace node from MPI packed buffer
 *            <buf>, starting from position <*pos>, where the total length
 *            of the buffer in bytes is <n>. 
 *
 * Returns:   <eslOK> on success. <*pos> is updated to the position of
 *            the next element in <buf> to unpack (if any). <*ret_cp9_tr>
 *            contains a newly allocated cp9trace, which the caller is 
 *            responsible for free'ing.
 *            
 * Throws:    <eslESYS> on an MPI call failure. <eslEMEM> on allocation failure.
 *            In either case, <*ret_cp9_tr> is <NULL>, and the state of <buf>
 *            and <*pos> is undefined and should be considered to be corrupted.
 */
int
cm_cp9trace_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, CP9trace_t **ret_cp9_tr)
{
  int status;
  CP9trace_t *cp9_tr = NULL;
  int cp9_tr_n;

  status = MPI_Unpack (buf, n, pos, &cp9_tr_n,        1, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  CP9AllocTrace(cp9_tr_n, &cp9_tr);
  cp9_tr->tlen = cp9_tr_n;

  status = MPI_Unpack (buf, n, pos, cp9_tr->statetype, cp9_tr_n, MPI_CHAR,  comm); if (status != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, cp9_tr->nodeidx,   cp9_tr_n, MPI_INT,   comm); if (status != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, cp9_tr->pos,       cp9_tr_n, MPI_INT,   comm); if (status != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  *ret_cp9_tr = cp9_tr;
  return eslOK;

 ERROR:
  if (cp9_tr != NULL) 
    CP9FreeTrace(cp9_tr);
  *ret_cp9_tr = NULL;
  return status;
}


/* Function:  cm_digitized_sq_MPIPackSize()
 * Synopsis:  Calculates number of bytes needed to pack a 
 *            the name and dsq of a digitized ESL_SQ.
 * Incept:    EPN, Mon Sep  3 19:22:05 2007
 *
 * Purpose:   Calculate an upper bound on the number of bytes
 *            that <cm_digitized_sq_MPIPack()> will need to pack
 *            a digitized sq in a packed MPI message in 
 *            communicator <comm>; return that number of bytes 
 *            in <*ret_n>. 
 *            
 *            Caller will generally use this result to determine how
 *            to allocate a buffer before starting to pack into it.
 *
 * Returns:   <eslOK> on success, and <*ret_n> contains the answer.
 *
 * Throws:    <eslESYS> if an MPI call fails, and <*ret_n> is set to 0. 
 *
 * Note:      The sizing calls here need to stay matched up with
 *            the calls in <cm_digitized_sq_MPIPack()>.
 */
int
cm_digitized_sq_MPIPackSize(const ESL_SQ *sq, MPI_Comm comm, int *ret_n)
{
  int status;
  int sz;
  int n = 0;

  /* contract check */
  if(sq->dsq == NULL) ESL_XEXCEPTION(eslESYS, "cm_digitized_sq_MPIPackSize, sq not digitized.");

  /* space for sq->n, sq->name, and sq->dsq only */
  status = MPI_Pack_size(1,     MPI_INT,  comm, &sz);              n += sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  status = esl_mpi_PackOptSize(sq->name, -1, MPI_CHAR, comm, &sz); n += sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  status = MPI_Pack_size((sq->n+2), MPI_BYTE, comm, &sz);          n += sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");

  *ret_n = n;
  return eslOK;

 ERROR:
  *ret_n = 0;
  return status;
}

/* Function:  cm_digitized_sq_MPIPack()
 * Synopsis:  Packs minimal part of a sq in digital form into MPI buffer.
 * Incept:    EPN, Mon Sep  3 19:28:30 2007
 *
 * Purpose:   Packs the essential info of a sq <sq> into an 
 *            MPI packed message buffer <buf> of length <n> bytes, 
 *            starting at byte position
 *            <*position>, for MPI communicator <comm>.
 *
 * Returns:   <eslOK> on success; <buf> now contains the
 *            packed <sq>, and <*position> is set to the byte
 *            immediately following the last byte of the results
 *            in <buf>. 
 *
 * Throws:    <eslESYS> if an MPI call fails; or <eslEMEM> if the
 *            buffer's length <n> is overflowed by trying to pack
 *            <rnode> into <buf>. In either case, the state of
 *            <buf> and <*position> is undefined, and both should
 *            be considered to be corrupted.
 *
 */
int
cm_digitized_sq_MPIPack(const ESL_SQ *sq, char *buf, int n, int *position, MPI_Comm comm)
{
  int status;

  /* contract check */
  if(sq->dsq == NULL) ESL_EXCEPTION(eslESYS, "cm_digitized_sq_MPIPackSize, sq not digitized.");

  ESL_DPRINTF2(("cm_digitized_sq_MPIPack(): ready.\n"));
  
  status = MPI_Pack((int *) &(sq->n),             1, MPI_INT, buf, n, position,  comm);  if (status != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  status = esl_mpi_PackOpt(sq->name, -1, MPI_CHAR, buf, n, position,  comm);             if (status != eslOK) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((ESL_DSQ *) sq->dsq, sq->n+2,MPI_BYTE, buf, n, position, comm);   if (status != 0)     ESL_EXCEPTION(eslESYS, "pack failed");

  ESL_DPRINTF2(("cm_digitized_sq_MPIPack(): done. Packed %d bytes into buffer of size %d\n", *position, n));

  if (*position > n) ESL_EXCEPTION(eslEMEM, "buffer overflow");
  return eslOK;
}

/* Function:  cm_digitized_sq_MPIUnpack()
 * Synopsis:  Unpacks essential part of a digitized sq from an MPI buffer.
 * Incept:    EPN, Mon Sep  3 19:32:00 2007
 *
 * Purpose:   Unpack a newly allocated sq from MPI packed buffer
 *            <buf>, starting from position <*pos>, where the total length
 *            of the buffer in bytes is <n>. 
 *
 * Returns:   <eslOK> on success. <*pos> is updated to the position of
 *            the next element in <buf> to unpack (if any). <*ret_tr>
 *            contains a newly allocated digitized sq, which the caller is 
 *            responsible for free'ing.
 *            
 * Throws:    <eslESYS> on an MPI call failure. <eslEMEM> on allocation failure.
 *            In either case, <*ret_tr> is <NULL>, and the state of <buf>
 *            and <*pos> is undefined and should be considered to be corrupted.
 */
int
cm_digitized_sq_MPIUnpack(const ESL_ALPHABET *abc, char *buf, int n, int *pos, MPI_Comm comm, ESL_SQ **ret_sq)
{
  int status;
  ESL_DSQ *dsq = NULL;
  char *name = NULL;
  int L;

  status = MPI_Unpack       (buf, n, pos, &L,                 1, MPI_INT,  comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  ESL_ALLOC(dsq, sizeof(ESL_DSQ) * (L+2));
  status = esl_mpi_UnpackOpt(buf, n, pos, (void**)&(name), NULL, MPI_CHAR, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack       (buf, n, pos, dsq,  L+2,  MPI_BYTE, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  dsq[0] = dsq[(L+1)] = eslDSQ_SENTINEL; /* overwrite */

  *ret_sq = esl_sq_CreateDigitalFrom(abc, name, dsq, L, NULL, NULL, NULL);
  free(dsq);  /* a copy was made */
  free(name); /* a copy was made */

  return eslOK;

 ERROR:
  if (dsq  != NULL) free(dsq);
  if (name != NULL) free(name);
  *ret_sq = NULL;
  return status;
}


/* Function:  exp_info_MPIPackSize()
 * Synopsis:  Calculates number of bytes needed to pack a 
 *            ExpInfo_t object. Follows 'Purpose' 
 *            of other *_MPIPackSize() functions above. 
 *
 * Incept:    EPN, Wed Dec 12 05:00:01 2007
 *
 * Returns:   <eslOK> on success, and <*ret_n> contains the answer.
 *
 * Throws:    <eslESYS> if an MPI call fails, and <*ret_n> is set to 0. 
 *
 * Note:      The sizing calls here need to stay matched up with
 *            the calls in <exp_info_MPIPack()>.
 */
int
exp_info_MPIPackSize(ExpInfo_t *exp, MPI_Comm comm, int *ret_n)
{
  int status;
  int sz;
  int n = 0;

  status = MPI_Pack_size(1, MPI_INT, comm, &sz);    n += 2*sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /* nrandhits, is_valid */
  status = MPI_Pack_size(1, MPI_LONG, comm, &sz);   n += 2*sz;  if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /* dbsize, cur_eff_dbsize */
  status = MPI_Pack_size(1, MPI_DOUBLE, comm, &sz); n += 4*sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /* mu_orig, mu_extrap, lambda, tailp */

  *ret_n = n;
  return eslOK;

 ERROR:
  *ret_n = 0;
  return status;
}

/* Function:  exp_info_MPIPack()
 * Synopsis:  Packs ExpInfo_t <exp> into MPI buffer.
 *            See 'Purpose','Returns' and 'Throws'
 *            of other *_MPIPack()'s for more info.
 *
 * Incept:    EPN, Wed Dec 12 05:03:40 2007
 */
int
exp_info_MPIPack(ExpInfo_t *exp, char *buf, int n, int *position, MPI_Comm comm)
{
  int status;

  ESL_DPRINTF2(("exp_info_MPIPack(): ready.\n"));
  
  status = MPI_Pack((long *)   &(exp->cur_eff_dbsize), 1, MPI_LONG,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((double *) &(exp->lambda),         1, MPI_DOUBLE, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((double *) &(exp->mu_extrap),      1, MPI_DOUBLE, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((double *) &(exp->mu_orig),        1, MPI_DOUBLE, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((long *)   &(exp->dbsize),         1, MPI_LONG,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *)    &(exp->nrandhits),      1, MPI_INT,    buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((double *) &(exp->tailp),          1, MPI_DOUBLE, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *)    &(exp->is_valid),       1, MPI_INT,    buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  ESL_DPRINTF2(("exp_info_results_MPIPack(): done. Packed %d bytes into buffer of size %d\n", *position, n));

  if (*position > n) ESL_EXCEPTION(eslEMEM, "buffer overflow");
  return eslOK;
}

/* Function:  exp_info_MPIUnpack()
 * Synopsis:  Unpacks ExpInfo_t <exp> from an MPI buffer.
 *            Follows 'Purpose', 'Returns', 'Throws' of other
 *            *_MPIUnpack() functions above.
 *
 * Incept:    EPN, Wed Dec 12 05:29:19 2007
 */
int
exp_info_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, ExpInfo_t **ret_exp)
{
  int status;
  ExpInfo_t *exp;

  ESL_ALLOC(exp, sizeof(ExpInfo_t));
  status = MPI_Unpack (buf, n, pos, &(exp->cur_eff_dbsize),  1, MPI_LONG,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &(exp->lambda),          1, MPI_DOUBLE, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &(exp->mu_extrap),       1, MPI_DOUBLE, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &(exp->mu_orig),         1, MPI_DOUBLE, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &(exp->dbsize),          1, MPI_LONG,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &(exp->nrandhits),       1, MPI_INT,    comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &(exp->tailp),           1, MPI_DOUBLE, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &(exp->is_valid),        1, MPI_INT,    comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  *ret_exp = exp;
  return eslOK;

 ERROR:
  if(exp != NULL) free(exp);
  *ret_exp = NULL;
  return status;
}

/* Function:  hmm_filter_info_MPIPackSize()
 * Synopsis:  Calculates number of bytes needed to pack a 
 *            HMMFilterInfo_t object. Follows 'Purpose' 
 *            of other *_MPIPackSize() functions above. 
 *
 * Incept:    EPN, Wed Dec 12 05:29:15 2007
 *
 * Returns:   <eslOK> on success, and <*ret_n> contains the answer.
 *
 * Throws:    <eslESYS> if an MPI call fails, and <*ret_n> is set to 0. 
 *
 * Note:      The sizing calls here need to stay matched up with
 *            the calls in <hmm_filter_info_MPIPack()>.
 */
int
hmm_filter_info_MPIPackSize(HMMFilterInfo_t *hfi, MPI_Comm comm, int *ret_n)
{
  int status;
  int sz;
  int n = 0;

  status = MPI_Pack_size(1, MPI_INT, comm, &sz);    n += 4*sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /*  is_valid, N, ncut, always_better_than_Smax */
  status = MPI_Pack_size(1, MPI_LONG,comm, &sz);    n +=  sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /*  dbsize */
  status = MPI_Pack_size(1, MPI_FLOAT, comm, &sz);  n += sz;   if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /*  F */
  status = MPI_Pack_size(1, MPI_INT, comm, &sz);    n += hfi->ncut*2*sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /* cm_E_cut, fwd_E_cut arrays */
  *ret_n = n;
  return eslOK;

 ERROR:
  *ret_n = 0;
  return status;
}

/* Function:  hmm_filter_info_MPIPack()
 * Synopsis:  Packs HMMFilterInfo_t <hfi> into MPI buffer.
 *            See 'Purpose','Returns' and 'Throws'
 *            of other *_MPIPack()'s for more info.
 *
 * Incept:    EPN, Wed Dec 12 05:29:10 2007
 */
int
hmm_filter_info_MPIPack(HMMFilterInfo_t *hfi, char *buf, int n, int *position, MPI_Comm comm)
{
  int status;

  ESL_DPRINTF2(("hmm_filter_info_MPIPack(): ready.\n"));
  
  status = MPI_Pack((int *)    &(hfi->is_valid),                1, MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *)    &(hfi->N),                       1, MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *)    &(hfi->dbsize),                  1, MPI_LONG,  buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *)    &(hfi->ncut),                    1, MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *)    &(hfi->always_better_than_Smax), 1, MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((float *)  &(hfi->F),                       1, MPI_FLOAT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((float *)  hfi->cm_E_cut,           hfi->ncut, MPI_FLOAT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((float *)  hfi->fwd_E_cut,          hfi->ncut, MPI_FLOAT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");

  if (*position > n) ESL_EXCEPTION(eslEMEM, "buffer overflow");
  return eslOK;
}

/* Function:  hmm_filter_info_MPIUnpack()
 * Synopsis:  Unpacks HMMFilterInfo_t <hfi> from an MPI buffer.
 *            Follows 'Purpose', 'Returns', 'Throws' of other
 *            *_MPIUnpack() functions above.
 *
 * Incept:    EPN, Wed Dec 12 05:29:05 2007
 */
int
hmm_filter_info_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, HMMFilterInfo_t **ret_hfi)
{
  int status;
  HMMFilterInfo_t *hfi;

  ESL_ALLOC(hfi, sizeof(HMMFilterInfo_t));
  status = MPI_Unpack (buf, n, pos, &(hfi->is_valid),                1, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &(hfi->N),                       1, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &(hfi->dbsize),                  1, MPI_LONG,  comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &(hfi->ncut),                    1, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &(hfi->always_better_than_Smax), 1, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &(hfi->F),                       1, MPI_FLOAT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  ESL_ALLOC(hfi->cm_E_cut,  sizeof(float) * hfi->ncut);
  ESL_ALLOC(hfi->fwd_E_cut, sizeof(float) * hfi->ncut);
  status = MPI_Unpack (buf, n, pos, hfi->cm_E_cut,           hfi->ncut, MPI_FLOAT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, hfi->fwd_E_cut,          hfi->ncut, MPI_FLOAT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  *ret_hfi = hfi;
  return eslOK;

 ERROR:
  if(hfi != NULL) FreeHMMFilterInfo(hfi);
  *ret_hfi = NULL;
  return status;
}

/* Function:  best_filter_info_MPIPackSize()
 * Synopsis:  Calculates number of bytes needed to pack a 
 *            BestFilterInfo_t object. Follows 'Purpose' 
 *            of other *_MPIPackSize() functions above. 
 *
 * Incept:    EPN, Wed Dec 12 05:29:15 2007
 *
 * Returns:   <eslOK> on success, and <*ret_n> contains the answer.
 *
 * Throws:    <eslESYS> if an MPI call fails, and <*ret_n> is set to 0. 
 *
 * Note:      The sizing calls here need to stay matched up with
 *            the calls in <best_filter_info_MPIPack()>.
 */
int
best_filter_info_MPIPackSize(BestFilterInfo_t *bf, MPI_Comm comm, int *ret_n)
{
  int status;
  int sz;
  int n = 0;

  status = MPI_Pack_size(1, MPI_INT, comm, &sz);    n += 6*sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /* cm_M, N, db_size, is_valid, ftype, np */
  status = MPI_Pack_size(1, MPI_FLOAT, comm, &sz);  n += 6*sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /* cm_eval, F, e_cutoff, full_cm_ncalcs, fil_ncalcs, fil_plus_surv_ncalcs */
  *ret_n = n;
  return eslOK;

 ERROR:
  *ret_n = 0;
  return status;
}

/* Function:  best_filter_info_MPIPack()
 * Synopsis:  Packs BestFilterInfo_t <bf> into MPI buffer.
 *            See 'Purpose','Returns' and 'Throws'
 *            of other *_MPIPack()'s for more info.
 *
 * Incept:    EPN, Wed Dec 12 05:29:10 2007
 */
int
best_filter_info_MPIPack(BestFilterInfo_t *bf, char *buf, int n, int *position, MPI_Comm comm)
{
  int status;

  ESL_DPRINTF2(("best_filter_info_MPIPack(): ready.\n"));
  
  status = MPI_Pack((int *) &(bf->cm_M),                  1, MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *) &(bf->N),                     1, MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *) &(bf->db_size),               1, MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *) &(bf->is_valid),              1, MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *) &(bf->ftype),                 1, MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");

  status = MPI_Pack((float *) &(bf->cm_eval),             1, MPI_FLOAT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((float *) &(bf->F),                   1, MPI_FLOAT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((float *) &(bf->e_cutoff),            1, MPI_FLOAT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((float *) &(bf->full_cm_ncalcs),      1, MPI_FLOAT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((float *) &(bf->fil_ncalcs),          1, MPI_FLOAT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((float *) &(bf->fil_plus_surv_ncalcs),1, MPI_FLOAT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");

  if (*position > n) ESL_EXCEPTION(eslEMEM, "buffer overflow");
  return eslOK;
}

/* Function:  best_filter_info_MPIUnpack()
 * Synopsis:  Unpacks BestFilterInfo_t <bf> from an MPI buffer.
 *            Follows 'Purpose', 'Returns', 'Throws' of other
 *            *_MPIUnpack() functions above.
 *
 * Incept:    EPN, Wed Dec 12 05:29:05 2007
 */
int
best_filter_info_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, BestFilterInfo_t **ret_bf)
{
  int status;
  BestFilterInfo_t *bf;

  ESL_ALLOC(bf, sizeof(BestFilterInfo_t));
  status = MPI_Unpack (buf, n, pos, &(bf->cm_M),                1, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &(bf->N),                   1, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &(bf->db_size),             1, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &(bf->is_valid),            1, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &(bf->ftype),               1, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  status = MPI_Unpack (buf, n, pos, &(bf->cm_eval),             1, MPI_FLOAT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &(bf->F),                   1, MPI_FLOAT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &(bf->e_cutoff),            1, MPI_FLOAT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &(bf->full_cm_ncalcs),      1, MPI_FLOAT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &(bf->fil_ncalcs),          1, MPI_FLOAT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &(bf->fil_plus_surv_ncalcs),1, MPI_FLOAT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  *ret_bf = bf;
  return eslOK;

 ERROR:
  if(bf != NULL) free(bf);
  *ret_bf = NULL;
  return status;
}



/* Function:  cmcalibrate_exp_results_MPIPackSize()
 * Synopsis:  Calculates number of bytes needed to pack a 
 *            results for CM scan for cmcalibrate.
 * Incept:    EPN, Thu Dec  6 16:56:27 2007
 *
 * Purpose:   Calculate an upper bound on the number of bytes
 *            that <cmcalibrate_exp_results_MPIPack()> will need 
 *            to pack it's results in a packed MPI message in 
 *            communicator <comm>; return that number of bytes 
 *            in <*ret_n>. 
 *            
 *            Caller will generally use this result to determine how
 *            to allocate a buffer before starting to pack into it.
 *
 * Returns:   <eslOK> on success, and <*ret_n> contains the answer.
 *
 * Throws:    <eslESYS> if an MPI call fails, and <*ret_n> is set to 0. 
 *
 * Note:      The sizing calls here need to stay matched up with
 *            the calls in <cmcalibrate_exp_results_MPIPack()>.
 */
int
cmcalibrate_exp_results_MPIPackSize(float *scA, int nseq, MPI_Comm comm, int *ret_n)
{
  int status;
  int sz;
  int n = 0;

  status = MPI_Pack_size(1, MPI_INT, comm, &sz);   n += sz;      if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  status = MPI_Pack_size(1, MPI_FLOAT, comm, &sz); n += nseq*sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");

  *ret_n = n;
  return eslOK;

 ERROR:
  *ret_n = 0;
  return status;
}

/* Function:  cmcalibrate_exp_results_MPIPack()
 * Synopsis:  Packs CM vscAA scores into MPI buffer.
 * Incept:    EPN, Thu Dec  6 16:56:31 2007
 *
 * Purpose:   Packs <vscAA> into an MPI packed message 
 *            buffer <buf> of length <n> bytes, 
 *            starting at byte position
 *            <*position>, for MPI communicator <comm>.
 *
 *            <scA> is an array, scA[0..i..nseq-1]
 *            holding the best score for a CM or HMM scan against sequence i.
 * 
 * Returns:   <eslOK> on success; <buf> now contains the
 *            packed <tr>, and <*position> is set to the byte
 *            immediately following the last byte of the results
 *            in <buf>. 
 *
 * Throws:    <eslESYS> if an MPI call fails; or <eslEMEM> if the
 *            buffer's length <n> is overflowed by trying to pack
 *            <rnode> into <buf>. In either case, the state of
 *            <buf> and <*position> is undefined, and both should
 *            be considered to be corrupted.
 *
 */
int
cmcalibrate_exp_results_MPIPack(float *scA, int nseq, char *buf, int n, int *position, MPI_Comm comm)
{
  int status;

  ESL_DPRINTF2(("cmcalibrate_exp_results_MPIPack(): ready.\n"));
  
  status = MPI_Pack((int *) &(nseq), 1, MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(scA,       nseq, MPI_FLOAT, buf, n, position,  comm);     if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  ESL_DPRINTF2(("cmcalibrate_exp_results_MPIPack(): done. Packed %d bytes into buffer of size %d\n", *position, n));

  if (*position > n) ESL_EXCEPTION(eslEMEM, "buffer overflow");
  return eslOK;
}

/* Function:  cmcalibrate_exp_results_MPIUnpack()
 * Synopsis:  Unpacks <scA> from an MPI buffer.
 * Incept:    EPN, Thu Dec  6 16:56:36 2007
 *
 * Purpose:   Unpack a newly allocated set of CM or HMM scores <scA> from MPI packed buffer
 *            <buf>, starting from position <*pos>, where the total length
 *            of the buffer in bytes is <n>. 
 *
 *            <scA> is an array, scA[0..i..nseq-1]
 *            holding the best score for a CM or HMM scan against sequence i.
 *
 * Returns:   <eslOK> on success. <*pos> is updated to the position of
 *            the next element in <buf> to unpack (if any). <*ret_tr>
 *            contains a newly allocated parsetree, which the caller is 
 *            responsible for free'ing.
 *            
 * Throws:    <eslESYS> on an MPI call failure. <eslEMEM> on allocation failure.
 *            In either case, <*ret_vscAA> is <NULL>, and the state of <buf>
 *            and <*pos> is undefined and should be considered to be corrupted.
 */
int
cmcalibrate_exp_results_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, float **ret_scA, int *ret_nseq)
{
  int status;
  float *scA;
  int nseq = 0;

  status = MPI_Unpack (buf, n, pos, &nseq,        1, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  ESL_ALLOC(scA, sizeof(float) * nseq);
  status = MPI_Unpack (buf, n, pos, scA, nseq, MPI_FLOAT,  comm);  if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  *ret_scA = scA;
  *ret_nseq   = nseq;
  return eslOK;

 ERROR:
  if(scA != NULL) free(scA);
  *ret_scA = NULL;
  *ret_nseq = 0;
  return status;
}

/* Function:  cmcalibrate_filter_results_MPIPackSize()
 * Synopsis:  Calculates number of bytes needed to pack 
 *            HMm filter results for cmcalibrate.
 *            Follows, 'Purpose', 'Returns', 'Throws' of
 *            the many other *_MPIPackSize() funcs above.
 *            
 * Incept:    EPN, Tue Jan  8 15:14:17 2008
 *           
 * Note:      The sizing calls here need to stay matched up with
 *            the calls in <cmcalibrate_filter_results_MPIPack()>.
 */
int
cmcalibrate_filter_results_MPIPackSize(int nseq, MPI_Comm comm, int *ret_n)
{
  int status;
  int sz;
  int n = 0;

  status = MPI_Pack_size(1, MPI_INT, comm, &sz);        n += sz;  if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /* for nseq */
  status = MPI_Pack_size(1, MPI_INT, comm, &sz);        n += sz;  if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /* for seq_offset */
  status = MPI_Pack_size(nseq, MPI_FLOAT, comm, &sz);   n += sz;  if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /* for cyk_scA */
  status = MPI_Pack_size(nseq, MPI_FLOAT, comm, &sz);   n += sz;  if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /* for ins_scA */
  status = MPI_Pack_size(nseq, MPI_FLOAT, comm, &sz);   n += sz;  if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /* for fwd_scA */

  *ret_n = n;
  return eslOK;

 ERROR:
  *ret_n = 0;
  return status;
}

/* Function:  cmcalibrate_filter_results_MPIPack()
 * Synopsis:  Packs cmcalibrate HMM filter results into MPI buffer.
 *            Follows, 'Purpose', 'Returns', 'Throws' of
 *            the many other *_MPIPack() funcs above.
 * 
 * Incept:    EPN, Wed Dec 12 16:36:02 2007
 *
 */
int
cmcalibrate_filter_results_MPIPack(float *cyk_scA, float *ins_scA, float *fwd_scA, int nseq, int seq_offset, char *buf, int n, int *position, MPI_Comm comm)
{
  int status;

  ESL_DPRINTF2(("cmcalibrate_filter_results_MPIPack(): ready.\n"));

  status = MPI_Pack((int *) &(nseq),       1, MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *) &(seq_offset), 1, MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(cyk_scA,            nseq, MPI_FLOAT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(ins_scA,            nseq, MPI_FLOAT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(fwd_scA,            nseq, MPI_FLOAT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");

  ESL_DPRINTF2(("cmcalibrate_filter_results_MPIPack(): done. Packed %d bytes into buffer of size %d\n", *position, n));

  if (*position > n) ESL_EXCEPTION(eslEMEM, "buffer overflow");
  return eslOK;

}

/* Function:  cmcalibrate_filter_results_MPIUnpack()
 * Synopsis:  Unpacks cmcalibrate HMM filter results from an MPI buffer.
 *            Follows, 'Purpose', 'Returns', 'Throws' of
 *            the many other *_MPIUnpack() funcs above.
 * Incept:    EPN, Wed Dec 12 16:38:15 2007
 *
 */
int
cmcalibrate_filter_results_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, float **ret_cyk_scA, float **ret_ins_scA, float **ret_fwd_scA, int *ret_nseq, int *ret_seq_offset)
{
  int status;
  float  *cyk_scA  = NULL;
  float  *ins_scA  = NULL;
  float  *fwd_scA  = NULL;
  int nseq = 0;
  int seq_offset = 0;

  status = MPI_Unpack (buf, n, pos, &nseq,        1, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &seq_offset,  1, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  ESL_ALLOC(cyk_scA, sizeof(float) * nseq);
  status = MPI_Unpack (buf, n, pos, cyk_scA, nseq, MPI_FLOAT,  comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  ESL_ALLOC(ins_scA, sizeof(float) * nseq);
  status = MPI_Unpack (buf, n, pos, ins_scA, nseq, MPI_FLOAT,  comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  ESL_ALLOC(fwd_scA, sizeof(float) * nseq);
  status = MPI_Unpack (buf, n, pos, fwd_scA, nseq, MPI_FLOAT,  comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  *ret_cyk_scA    = cyk_scA;
  *ret_ins_scA    = ins_scA;
  *ret_fwd_scA    = fwd_scA;
  *ret_nseq       = nseq;
  *ret_seq_offset = seq_offset;
  return eslOK;

  ESL_DPRINTF1(("cmcalibrate_filter_results_MPIUnpack() done.\n"));

 ERROR:
  if(cyk_scA != NULL) free(cyk_scA);
  if(ins_scA != NULL) free(ins_scA);
  if(fwd_scA != NULL) free(fwd_scA);
  *ret_nseq = 0;
  *ret_seq_offset = 0;
  return status;
}

/* Function:  comlog_MPIPackSize()
 * Synopsis:  Calculates number of bytes needed to pack a 
 *            ComLog_t object. Follows 'Purpose' 
 *            of other *_MPIPackSize() functions above. 
 *
 * Incept:    EPN, Mon Dec 31 14:31:04 2007
 *
 * Returns:   <eslOK> on success, and <*ret_n> contains the answer.
 *
 * Throws:    <eslESYS> if an MPI call fails, and <*ret_n> is set to 0. 
 *
 * Note:      The sizing calls here need to stay matched up with
 *            the calls in <comlog_MPIPack()>.
 */
int
comlog_MPIPackSize(ComLog_t *comlog, MPI_Comm comm, int *ret_n)
{
  int status;
  int sz;
  int n = 0;

  status = esl_mpi_PackOptSize(comlog->bcom,  -1, MPI_CHAR, comm, &sz); n += sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  status = esl_mpi_PackOptSize(comlog->bdate, -1, MPI_CHAR, comm, &sz); n += sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  status = esl_mpi_PackOptSize(comlog->ccom,  -1, MPI_CHAR, comm, &sz); n += sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  status = esl_mpi_PackOptSize(comlog->cdate, -1, MPI_CHAR, comm, &sz); n += sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");

  *ret_n = n;
  return eslOK;

 ERROR:
  return status;
}

/* Function:  comlog_MPIPack()
 * Synopsis:  Packs ComLog_t <comlog> into MPI buffer.
 *            See 'Purpose','Returns' and 'Throws'
 *            of other *_MPIPack()'s for more info.
 *
 * Incept:    EPN, Wed Dec 12 05:03:40 2007
 */
int
comlog_MPIPack(ComLog_t *comlog, char *buf, int n, int *position, MPI_Comm comm)
{
  int status;

  ESL_DPRINTF2(("comlog_MPIPack(): ready.\n"));
  
  status = esl_mpi_PackOpt(comlog->bcom,   -1, MPI_CHAR, buf, n, position,  comm); if (status != eslOK) ESL_EXCEPTION(eslESYS, "pack failed");
  status = esl_mpi_PackOpt(comlog->bdate,  -1, MPI_CHAR, buf, n, position,  comm); if (status != eslOK) ESL_EXCEPTION(eslESYS, "pack failed");
  status = esl_mpi_PackOpt(comlog->ccom,   -1, MPI_CHAR, buf, n, position,  comm); if (status != eslOK) ESL_EXCEPTION(eslESYS, "pack failed");
  status = esl_mpi_PackOpt(comlog->cdate,  -1, MPI_CHAR, buf, n, position,  comm); if (status != eslOK) ESL_EXCEPTION(eslESYS, "pack failed");
  ESL_DPRINTF2(("comlog_results_MPIPack(): done. Packed %d bytes into buffer of size %d\n", *position, n));

  if (*position > n) ESL_EXCEPTION(eslEMEM, "buffer overflow");
  return eslOK;
}

/* Function:  comlog_MPIUnpack()
 * Synopsis:  Unpacks ComLog_t<comlog> from an MPI buffer.
 *            Follows 'Purpose', 'Returns', 'Throws' of other
 *            *_MPIUnpack() functions above.
 *
 * Incept:    EPN, Wed Dec 12 05:29:19 2007
 */
int
comlog_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, ComLog_t **ret_comlog)
{
  int status;
  ComLog_t *comlog;

  comlog = CreateComLog();
  status = esl_mpi_UnpackOpt(buf, n, pos, (void**)&(comlog->bcom),  NULL, MPI_CHAR, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = esl_mpi_UnpackOpt(buf, n, pos, (void**)&(comlog->bdate), NULL, MPI_CHAR, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = esl_mpi_UnpackOpt(buf, n, pos, (void**)&(comlog->ccom),  NULL, MPI_CHAR, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = esl_mpi_UnpackOpt(buf, n, pos, (void**)&(comlog->cdate), NULL, MPI_CHAR, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  *ret_comlog = comlog;
  return eslOK;

 ERROR:
  if(comlog != NULL) free(comlog);
  *ret_comlog = NULL;
  return status;
}

/*****************************************************************
 * Communicating CM_PIPELINE
 *****************************************************************/

/* Function:  cm_pipeline_MPISend()
 * Synopsis:  Send pipeline data as an MPI message.
 * Incept:    EPN, Thu Jun  2 14:36:01 2011
 *            MSF, Wed Sep 09 09:09:47 2009 [Janelia] (p7_pipeline_MPISend())
 *
 * Purpose:   Sends pipeline statistics <pli> to MPI process <dest>
 *            (where <dest> ranges from 0..<nproc-1>), with MPI tag
 *            <tag> for MPI communicator <comm>.
 *            
 *            In order to minimize alloc/free cycles in this routine,
 *            caller passes a pointer to a working buffer <*buf> of
 *            size <*nalloc> characters. If necessary (i.e. if <pli> is
 *            too big to fit), <*buf> will be reallocated and <*nalloc>
 *            increased to the new size. As a special case, if <*buf>
 *            is <NULL> and <*nalloc> is 0, the buffer will be
 *            allocated appropriately, but the caller is still
 *            responsible for free'ing it.
 *            
 *            If <pli> is NULL, the pipeline statistics are initialized
 *            to zeros.
 *            
 * Returns:   <eslOK> on success.
 */
int
cm_pipeline_MPISend(CM_PIPELINE *pli, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc)
{
  int   status;
  int   sz, n, pos;

  CM_PIPELINE bogus;

  /* This will look wasteful, but the MPI spec doesn't guarantee that 
   * MPI_Pack_size(x, ...) + MPI_Pack_size(x, ... ) == MPI_Pack_size(2x, ...).
   * Indeed there are some hints in the spec that that's *not* true.
   * So we assume we must match our Pack_size calls exactly to our Pack calls.
   */
  n = 0;
  if (MPI_Pack_size(1, MPI_LONG_INT,      comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  n += sz; /* mode */
  if (MPI_Pack_size(1, MPI_LONG_INT,      comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  n += sz; /* Z_setby */
  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  n += sz; /* nmodels */
  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  n += sz; /* nseqs */
  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  n += sz; /* nres */
  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  n += sz; /* nnodes */
  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  n += sz; /* n_past_msv */
  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  n += sz; /* n_past_vit */
  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  n += sz; /* n_past_fwd */
  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  n += sz; /* n_past_gfwd */
  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  n += sz; /* n_past_edef */
  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  n += sz; /* n_past_cyk */
  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  n += sz; /* n_past_ins */
  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  n += sz; /* n_past_msvbias */
  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  n += sz; /* n_past_vitbias */
  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  n += sz; /* n_past_fwdbias */
  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  n += sz; /* n_past_gfwdbias */
  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  n += sz; /* n_past_edefbias */
  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  n += sz; /* pos_past_msv */
  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  n += sz; /* pos_past_vit */
  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  n += sz; /* pos_past_fwd */
  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  n += sz; /* pos_past_gfwd */
  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  n += sz; /* pos_past_edef */
  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  n += sz; /* pos_past_cyk */
  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  n += sz; /* pos_past_ins */
  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  n += sz; /* pos_past_msvbias */
  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  n += sz; /* pos_past_vitbias */
  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  n += sz; /* pos_past_fwdbias */
  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  n += sz; /* pos_past_gfwdbias */
  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  n += sz; /* pos_past_edefbias */

  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  n += sz; /* n_overflow_fcyk  */
  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  n += sz; /* n_overflow_final */
  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  n += sz; /* n_aln_hboa       */
  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  n += sz; /* n_aln_hbcyk      */
  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  n += sz; /* n_aln_dccyk      */

  if (MPI_Pack_size(1, MPI_DOUBLE,        comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");  n += sz; /* Z */
  
  /* Make sure the buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) {
    void *tmp;
    ESL_RALLOC(*buf, tmp, sizeof(char) * n);
    *nalloc = n; 
  }

  /* if no pipeline was defined, return zeros for the stats */
  if (pli == NULL) 
    {
      bogus.mode              = CM_SEARCH_SEQS;    /* that's 0. (some compilers complain if you set 0 directly. */
      bogus.Z_setby           = CM_ZSETBY_SSIINFO; /* ditto. */
      bogus.nmodels           = 0;
      bogus.nseqs             = 0;
      bogus.nres              = 0;
      bogus.nnodes            = 0;
      bogus.n_past_msv        = 0;
      bogus.n_past_vit        = 0;
      bogus.n_past_fwd        = 0;
      bogus.n_past_gfwd       = 0;
      bogus.n_past_edef       = 0;
      bogus.n_past_cyk        = 0;
      bogus.n_past_ins        = 0;
      bogus.n_past_msvbias    = 0;
      bogus.n_past_vitbias    = 0;
      bogus.n_past_fwdbias    = 0;
      bogus.n_past_gfwdbias   = 0;
      bogus.n_past_edefbias   = 0;
      bogus.pos_past_msv      = 0;
      bogus.pos_past_vit      = 0;
      bogus.pos_past_fwd      = 0;
      bogus.pos_past_gfwd     = 0;
      bogus.pos_past_edef     = 0;
      bogus.pos_past_cyk      = 0;
      bogus.pos_past_ins      = 0;
      bogus.pos_past_msvbias  = 0;
      bogus.pos_past_vitbias  = 0;
      bogus.pos_past_fwdbias  = 0;
      bogus.pos_past_gfwdbias = 0;
      bogus.pos_past_edefbias = 0;
      bogus.n_overflow_fcyk   = 0;
      bogus.n_overflow_final  = 0;
      bogus.n_aln_hboa        = 0;
      bogus.n_aln_hbcyk       = 0;
      bogus.n_aln_dccyk       = 0;
      bogus.Z                 = 0.0;
      pli = &bogus;
   } 

  /* Pack the pipeline into the buffer */
  pos = 0;
  if (MPI_Pack(&pli->mode,            1, MPI_LONG_INT,      *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&pli->Z_setby,         1, MPI_LONG_INT,      *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&pli->nmodels,         1, MPI_LONG_LONG_INT, *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&pli->nseqs,           1, MPI_LONG_LONG_INT, *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&pli->nres,            1, MPI_LONG_LONG_INT, *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&pli->nnodes,          1, MPI_LONG_LONG_INT, *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&pli->n_past_msv,      1, MPI_LONG_LONG_INT, *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&pli->n_past_vit,      1, MPI_LONG_LONG_INT, *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&pli->n_past_fwd,      1, MPI_LONG_LONG_INT, *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&pli->n_past_gfwd,     1, MPI_LONG_LONG_INT, *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&pli->n_past_edef,     1, MPI_LONG_LONG_INT, *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&pli->n_past_cyk,      1, MPI_LONG_LONG_INT, *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&pli->n_past_ins,      1, MPI_LONG_LONG_INT, *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&pli->n_past_msvbias,  1, MPI_LONG_LONG_INT, *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&pli->n_past_vitbias,  1, MPI_LONG_LONG_INT, *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&pli->n_past_fwdbias,  1, MPI_LONG_LONG_INT, *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&pli->n_past_gfwdbias, 1, MPI_LONG_LONG_INT, *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&pli->n_past_edefbias, 1, MPI_LONG_LONG_INT, *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 

  if (MPI_Pack(&pli->pos_past_msv,      1, MPI_LONG_LONG_INT, *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&pli->pos_past_vit,      1, MPI_LONG_LONG_INT, *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&pli->pos_past_fwd,      1, MPI_LONG_LONG_INT, *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&pli->pos_past_gfwd,     1, MPI_LONG_LONG_INT, *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&pli->pos_past_edef,     1, MPI_LONG_LONG_INT, *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&pli->pos_past_cyk,      1, MPI_LONG_LONG_INT, *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&pli->pos_past_ins,      1, MPI_LONG_LONG_INT, *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&pli->pos_past_msvbias,  1, MPI_LONG_LONG_INT, *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&pli->pos_past_vitbias,  1, MPI_LONG_LONG_INT, *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&pli->pos_past_fwdbias,  1, MPI_LONG_LONG_INT, *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&pli->pos_past_gfwdbias, 1, MPI_LONG_LONG_INT, *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&pli->pos_past_edefbias, 1, MPI_LONG_LONG_INT, *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 

  if (MPI_Pack(&pli->n_overflow_fcyk,   1, MPI_LONG_LONG_INT, *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&pli->n_overflow_final,  1, MPI_LONG_LONG_INT, *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&pli->n_aln_hboa,        1, MPI_LONG_LONG_INT, *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&pli->n_aln_hbcyk,       1, MPI_LONG_LONG_INT, *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&pli->n_aln_dccyk,       1, MPI_LONG_LONG_INT, *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 

  if (MPI_Pack(&pli->Z,                 1, MPI_DOUBLE,        *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 

  /* Send the packed pipeline to destination  */
  MPI_Send(*buf, n, MPI_PACKED, dest, tag, comm);
  return eslOK;
  
 ERROR:
  return status;
}


/* Function:  cm_pipeline_MPIRecv()
 * Synopsis:  Receive pipeline data as an MPI message.
 * Incept:    EPN, Thu Jun  2 14:41:46 2011
 *            MSF, Wed Sep 09 09:09:47 2009 [Janelia] (p7_pipeline_MPIRecv()
 *
 * Purpose:   Receive a pipeline from <source> (where <source> is usually
 *            process 0, the master) with tag <tag> from communicator <comm>,
 *            and return it in <*ret_pli>. 
 *            
 *            To minimize alloc/free cycles in this routine, caller
 *            passes a pointer to a buffer <*buf> of size <*nalloc>
 *            characters. These are passed by reference because if
 *            necessary, <buf> will be reallocated and <nalloc>
 *            increased to the new size. As a special case, if <buf>
 *            is <NULL> and <nalloc> is 0, the buffer will be
 *            allocated appropriately, but the caller is still
 *            responsible for free'ing it.
 *
 * Returns:   <eslOK> on success. <*ret_pli> contains the new pipeline;
 *            it is allocated here, and the caller is responsible for
 *            free'ing it.  <*buf> may have been reallocated to a
 *            larger size, and <*nalloc> may have been increased.
 *            
 */
int
cm_pipeline_MPIRecv(int source, int tag, MPI_Comm comm, char **buf, int *nalloc, ESL_GETOPTS *go, CM_PIPELINE **ret_pli)
{
  int          status;
  CM_PIPELINE *pli    = NULL;
  int          n;
  int          pos;
  MPI_Status   mpistatus;

  /* Probe first, because we need to know if our buffer is big enough.
   */
  MPI_Probe(source, tag, comm, &mpistatus);
  MPI_Get_count(&mpistatus, MPI_PACKED, &n);

  /* Make sure the buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) {
    void *tmp;
    ESL_RALLOC(*buf, tmp, sizeof(char) * n); 
    *nalloc = n; 
  }

  /* Receive the packed pipeline */
  MPI_Recv(*buf, n, MPI_PACKED, source, tag, comm, &mpistatus);

  /* Unpack it - watching out for the EOD signal of M = -1. */
  pos = 0;                                /* irrelevant, these are overwritten below */
  if ((pli = cm_pipeline_Create(go, 1, 1, 0, CM_ZSETBY_SSIINFO, CM_SEARCH_SEQS)) == NULL) { status = eslEMEM; goto ERROR; } /* mode will be immediately overwritten */
  if (MPI_Unpack(*buf, n, &pos, &(pli->mode),            1, MPI_LONG_INT,      comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(*buf, n, &pos, &(pli->Z_setby),         1, MPI_LONG_INT,      comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(*buf, n, &pos, &(pli->nmodels),         1, MPI_LONG_LONG_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(*buf, n, &pos, &(pli->nseqs),           1, MPI_LONG_LONG_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(*buf, n, &pos, &(pli->nres),            1, MPI_LONG_LONG_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(*buf, n, &pos, &(pli->nnodes),          1, MPI_LONG_LONG_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(*buf, n, &pos, &(pli->n_past_msv),      1, MPI_LONG_LONG_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(*buf, n, &pos, &(pli->n_past_vit),      1, MPI_LONG_LONG_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(*buf, n, &pos, &(pli->n_past_fwd),      1, MPI_LONG_LONG_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(*buf, n, &pos, &(pli->n_past_gfwd),     1, MPI_LONG_LONG_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(*buf, n, &pos, &(pli->n_past_edef),     1, MPI_LONG_LONG_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(*buf, n, &pos, &(pli->n_past_cyk),      1, MPI_LONG_LONG_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(*buf, n, &pos, &(pli->n_past_ins),      1, MPI_LONG_LONG_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(*buf, n, &pos, &(pli->n_past_msvbias),  1, MPI_LONG_LONG_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(*buf, n, &pos, &(pli->n_past_vitbias),  1, MPI_LONG_LONG_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(*buf, n, &pos, &(pli->n_past_fwdbias),  1, MPI_LONG_LONG_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(*buf, n, &pos, &(pli->n_past_gfwdbias), 1, MPI_LONG_LONG_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(*buf, n, &pos, &(pli->n_past_edefbias), 1, MPI_LONG_LONG_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed"); 

  if (MPI_Unpack(*buf, n, &pos, &(pli->pos_past_msv),      1, MPI_LONG_LONG_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(*buf, n, &pos, &(pli->pos_past_vit),      1, MPI_LONG_LONG_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(*buf, n, &pos, &(pli->pos_past_fwd),      1, MPI_LONG_LONG_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(*buf, n, &pos, &(pli->pos_past_gfwd),     1, MPI_LONG_LONG_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(*buf, n, &pos, &(pli->pos_past_edef),     1, MPI_LONG_LONG_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(*buf, n, &pos, &(pli->pos_past_cyk),      1, MPI_LONG_LONG_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(*buf, n, &pos, &(pli->pos_past_ins),      1, MPI_LONG_LONG_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(*buf, n, &pos, &(pli->pos_past_msvbias),  1, MPI_LONG_LONG_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(*buf, n, &pos, &(pli->pos_past_vitbias),  1, MPI_LONG_LONG_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(*buf, n, &pos, &(pli->pos_past_fwdbias),  1, MPI_LONG_LONG_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(*buf, n, &pos, &(pli->pos_past_gfwdbias), 1, MPI_LONG_LONG_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(*buf, n, &pos, &(pli->pos_past_edefbias), 1, MPI_LONG_LONG_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed"); 

  if (MPI_Unpack(*buf, n, &pos, &(pli->Z),           1, MPI_DOUBLE,        comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed"); 

  *ret_pli = pli;
  return eslOK;

 ERROR:
  if (pli != NULL) cm_pipeline_Destroy(pli, NULL);
  *ret_pli = NULL;
  return status;
}

/*****************************************************************
 * Communicating CM_TOPHITS
 *****************************************************************/

/* Function:  cm_tophits_MPISend()
 * Synopsis:  Send the TOPHITS as an MPI work unit.
 * Incept:    EPN, Thu Jun  2 08:45:05 2011
 *            MSF, Mon Sep  21 22:22:00 2009 [Janelia] (p7_tophits_MPISend()
 *
 * Purpose:   Sends the TOPHITS <th> as a work unit to MPI process
 *            <dest> (where <dest> ranges from 0..<nproc-1>), tagged
 *            with MPI tag <tag>, for MPI communicator <comm>, as 
 *            the sole workunit or result. 
 *            
 *            After the TOPHITS <th> information has been sent, send
 *            the each hit as an indepentant message.
 *            
 * Returns:   <eslOK> on success; <*buf> may have been reallocated and
 *            <*nalloc> may have been increased.
 * 
 * Throws:    <eslESYS> if an MPI call fails; <eslEMEM> if a malloc/realloc
 *            fails. In either case, <*buf> and <*nalloc> remain valid and useful
 *            memory (though the contents of <*buf> are undefined). 
 */
int
cm_tophits_MPISend(CM_TOPHITS *th, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc)
{
  int   status;
  int   sz, n, pos;
  int   i;

  CM_HIT    *hit = NULL;

  sz = 0;
  n  = 0;
  
  /* calculate the buffer size needed to hold the largest hit */
  hit = th->unsrt;
  for (i = 0; i < th->N; i++) {
    if ((status = cm_hit_MPIPackSize(hit, comm, &sz)) != eslOK) goto ERROR;
    n = ESL_MAX(n, sz);
    hit++;
  }
  if (MPI_Pack_size(3, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n = ESL_MAX(n, sz);

  /* Make sure the buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) {
    void *tmp;
    ESL_RALLOC(*buf, tmp, sizeof(char) * n);
    *nalloc = n; 
  }

  pos = 0;
  if (MPI_Pack(&th->N,         1, MPI_LONG_LONG_INT, *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&th->nreported, 1, MPI_LONG_LONG_INT, *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&th->nincluded, 1, MPI_LONG_LONG_INT, *buf, n, &pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 

  /* Send the packed tophits information */
  if (MPI_Send(*buf, n, MPI_PACKED, dest, tag, comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi send failed");
  if (th->N == 0) return eslOK;

  /* loop through the hit list sending to dest */
  hit = th->unsrt;
  for (i = 0; i < th->N; i++) {
    if ((status = cm_hit_MPISend(hit, dest, tag, comm, buf, nalloc)) != eslOK) goto ERROR;
    hit++;
  }

  return eslOK;

 ERROR:
  return status;
}

/* Function:  cm_tophits_MPIRecv()
 * Synopsis:  Receives an TOPHITS as a work unit from an MPI sender.
 * Incept:    EPN, Thu Jun  2 14:07:19 2011
 *            MSF, Mon Sep  21 22:22:00 2009 [Janelia] (p7_tophits_MPIRecv())
 *
 * Purpose:   Sends the TOPHITS <th> as a work unit to MPI process
 *            <dest> (where <dest> ranges from 0..<nproc-1>), tagged
 *            with MPI tag <tag>, for MPI communicator <comm>, as 
 *            the sole workunit or result. 
 *            
 *            After the TOPHITS <th> information has been sent, send
 *            the each hit as an indepentant message.
 *            
 * Returns:   <eslOK> on success; <*buf> may have been reallocated and
 *            <*nalloc> may have been increased.
 * 
 * Throws:    <eslESYS> if an MPI call fails; <eslEMEM> if a malloc/realloc
 *            fails. In either case, <*buf> and <*nalloc> remain valid and useful
 *            memory (though the contents of <*buf> are undefined). 
 */
int
cm_tophits_MPIRecv(int source, int tag, MPI_Comm comm, char **buf, int *nalloc, CM_TOPHITS **ret_th)
{
  int         n;
  int         status;
  int         pos;
  CM_TOPHITS *th    = NULL;
  CM_HIT     *hit   = NULL;
  MPI_Status  mpistatus;

  uint64_t    nhits;
  uint64_t    inx;

  /* Probe first, because we need to know if our buffer is big enough.
   */
  MPI_Probe(source, tag, comm, &mpistatus);
  MPI_Get_count(&mpistatus, MPI_PACKED, &n);

  /* make sure we are getting the tag we expect and from whom we expect if from */
  if (tag    != MPI_ANY_TAG    && mpistatus.MPI_TAG    != tag) {
    status = eslFAIL;
    goto ERROR;
  }
  if (source != MPI_ANY_SOURCE && mpistatus.MPI_SOURCE != source) {
    status = eslFAIL;
    goto ERROR;
  }

  /* set the source and tag */
  tag = mpistatus.MPI_TAG;
  source = mpistatus.MPI_SOURCE;

  /* Make sure the buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) {
    void *tmp;
    ESL_RALLOC(*buf, tmp, sizeof(char) * n); 
    *nalloc = n; 
  }

  /* Receive the packed top hits */
  MPI_Recv(*buf, n, MPI_PACKED, source, tag, comm, &mpistatus);

  /* Unpack it - watching out for the EOD signal of M = -1. */
  pos = 0;
  if ((th = cm_tophits_Create()) == NULL) { status = eslEMEM; goto ERROR; }
  if (MPI_Unpack(*buf, n, &pos, &nhits,         1, MPI_LONG_LONG_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(*buf, n, &pos, &th->nreported, 1, MPI_LONG_LONG_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(*buf, n, &pos, &th->nincluded, 1, MPI_LONG_LONG_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed");

  /* loop through all of the hits sent */
  for (inx = 0; inx < nhits; inx++) {
    if ((status = cm_tophits_CreateNextHit(th, &hit))                  != eslOK) goto ERROR;
    if ((status = cm_hit_MPIRecv(source, tag, comm, buf, nalloc, hit)) != eslOK) goto ERROR;
  }
  
  *ret_th = th;
  return eslOK;

 ERROR:
  if (th  != NULL) cm_tophits_Destroy(th);
  *ret_th = NULL;
  return status;
}


/* Function:  cm_hit_MPISend()
 */
int
cm_hit_MPISend(CM_HIT *hit, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc)
{
  int   status;
  int   pos;
  int   n = *nalloc;

  /* Pack the HIT into the buffer */
  pos  = 0;
  if ((status = cm_hit_MPIPack(hit, *buf, n, &pos, comm)) != eslOK) goto ERROR;

  /* Send the packed HIT to the destination. */
  if (MPI_Send(*buf, n, MPI_PACKED, dest, tag, comm) != 0)  ESL_EXCEPTION(eslESYS, "mpi send failed");

  return eslOK;

 ERROR:
  return status;
}

/* Function:  cm_hit_MPIRecv()
 */
int
cm_hit_MPIRecv(int source, int tag, MPI_Comm comm, char **buf, int *nalloc, CM_HIT *hit)
{
  int         n;
  int         status;
  int         pos;
  MPI_Status  mpistatus;

  /* Probe first, because we need to know if our buffer is big enough.
   */
  MPI_Probe(source, tag, comm, &mpistatus);
  MPI_Get_count(&mpistatus, MPI_PACKED, &n);

  /* make sure we are getting the tag we expect and from whom we expect if from */
  if (tag    != MPI_ANY_TAG    && mpistatus.MPI_TAG    != tag) {
    status = eslFAIL;
    goto ERROR;
  }
  if (source != MPI_ANY_SOURCE && mpistatus.MPI_SOURCE != source) {
    status = eslFAIL;
    goto ERROR;
  }

  /* set the source and tag */
  tag = mpistatus.MPI_TAG;
  source = mpistatus.MPI_SOURCE;

  /* Make sure the buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) {
    void *tmp;
    ESL_RALLOC(*buf, tmp, sizeof(char) * n); 
    *nalloc = n; 
  }

  /* Receive the packed top hits */
  MPI_Recv(*buf, n, MPI_PACKED, source, tag, comm, &mpistatus);

  /* Unpack it - watching out for the EOD signal of M = -1. */
  pos = 0;
  if ((status = cm_hit_MPIUnpack(*buf, n, &pos, comm, hit)) != eslOK) goto ERROR;
  
  return eslOK;

 ERROR:
  return status;
}

/* Function:  cm_hit_MPIPackSize()
 * Synopsis:  Calculates size needed to pack a HIT.
 * Incept:    MSF, Mon Sep  21 22:22:00 2009 [Janelia]
 *
 * Purpose:   Calculate an upper bound on the number of bytes
 *            that <cm_hit_MPIPack()> will need to pack an CM_HIT
 *            <hit> in a packed MPI message for MPI communicator
 *            <comm>; return that number of bytes in <*ret_n>.
 *
 * Returns:   <eslOK> on success, and <*ret_n> contains the answer.
 *
 * Throws:    <eslESYS> if an MPI call fails, and <*ret_n> is 0.
 */
int
cm_hit_MPIPackSize(CM_HIT *hit, MPI_Comm comm, int *ret_n)
{
  int   status;
  int   n = 0;
  int   sz;

  CM_ALIDISPLAY *ad = hit->ad;

  /* CM_HIT data */
  if (MPI_Pack_size(3,            MPI_LONG,   comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed"); n += sz;  /* start, stop, seq_idx    */
  if (MPI_Pack_size(1,            MPI_INT,    comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed"); n += sz;  /* in_rc                   */
  if (MPI_Pack_size(1,            MPI_FLOAT,  comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed"); n += sz;  /* score                   */
  if (MPI_Pack_size(2,            MPI_DOUBLE, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed"); n += sz;  /* pvalue, evalue          */
  if (MPI_Pack_size(1,            MPI_INT,    comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed"); n += sz;  /* flags                   */

  if ((status = esl_mpi_PackOptSize(hit->name, -1, MPI_CHAR, comm, &sz)) != eslOK) goto ERROR;  n += sz;
  if ((status = esl_mpi_PackOptSize(hit->acc,  -1, MPI_CHAR, comm, &sz)) != eslOK) goto ERROR;  n += sz; 
  if ((status = esl_mpi_PackOptSize(hit->desc, -1, MPI_CHAR, comm, &sz)) != eslOK) goto ERROR;  n += sz; 

  /* CM_ALIDISPLAY data */
  if (MPI_Pack_size(13,          MPI_INT,     comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed"); n += sz;  /* offset info             */
  if (MPI_Pack_size(4,           MPI_INT,     comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed"); n += sz;  /* N, cfrom, cto, clen     */
  if (MPI_Pack_size(3,           MPI_LONG,    comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed"); n += sz;  /* sequence info           */
  if (MPI_Pack_size(2,           MPI_INT,     comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed"); n += sz;  /* used_optacc, used_hbands*/
  if (MPI_Pack_size(2,           MPI_FLOAT,   comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed"); n += sz;  /* aln_sc, matrix_Mb       */
  if (MPI_Pack_size(1,           MPI_DOUBLE,  comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed"); n += sz;  /* elapsed_secs            */
  if (MPI_Pack_size(1,           MPI_INT,     comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed"); n += sz;  /* memsize                 */
  if (MPI_Pack_size(ad->memsize, MPI_CHAR,    comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed"); n += sz;  /* mem                     */

  *ret_n = n;
  return eslOK;

 ERROR:
  *ret_n = 0;
  return status;

}

/* Function:  cm_hit_MPIPack()
 * Synopsis:  Packs the HIT into MPI buffer.
 * Incept:    MSF, Mon Sep  21 22:22:00 2009 [Janelia]
 *
 * Purpose:   Packs HIT <hit> into an MPI packed message buffer <buf>
 *            of length <n> bytes, starting at byte position <*position>,
 *            for MPI communicator <comm>.
 *            
 *            The caller must know that <buf>'s allocation of <n>
 *            bytes is large enough to append the packed HIT at
 *            position <*pos>. This typically requires a call to
 *            <cm_hit_MPIPackSize()> first, and reallocation if
 *            needed.
 *            
 * Returns:   <eslOK> on success; <buf> now contains the
 *            packed <hit>, and <*position> is set to the byte
 *            immediately following the last byte of the HIT
 *            in <buf>. 
 *
 * Throws:    <eslESYS> if an MPI call fails; or <eslEMEM> if the
 *            buffer's length <n> was overflowed in trying to pack
 *            <msa> into <buf>. In either case, the state of
 *            <buf> and <*position> is undefined, and both should
 *            be considered to be corrupted.
 */
int
cm_hit_MPIPack(CM_HIT *hit, char *buf, int n, int *pos, MPI_Comm comm)
{
  int             status;
  int             offset;

  CM_ALIDISPLAY *ad = hit->ad;

  if (MPI_Pack(&hit->start,          1, MPI_LONG,     buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&hit->stop,           1, MPI_LONG,     buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&hit->in_rc,          1, MPI_INT,      buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&hit->seq_idx,        1, MPI_LONG,     buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&hit->score,          1, MPI_FLOAT,    buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&hit->pvalue,         1, MPI_DOUBLE,   buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&hit->evalue,         1, MPI_DOUBLE,   buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&hit->flags,          1, MPI_INT,      buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed");

  if ((status = esl_mpi_PackOpt(hit->name,        -1,      MPI_CHAR,  buf, n, pos, comm)) != eslOK) return status;
  if ((status = esl_mpi_PackOpt(hit->acc,         -1,      MPI_CHAR,  buf, n, pos, comm)) != eslOK) return status; 
  if ((status = esl_mpi_PackOpt(hit->desc,        -1,      MPI_CHAR,  buf, n, pos, comm)) != eslOK) return status; 

  offset = (ad->rfline  == NULL)  ? -1 : ad->rfline - ad->mem;
  if (MPI_Pack(&offset,              1, MPI_INT,      buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed");
  offset = (ad->nline  == NULL)   ? -1 : ad->nline - ad->mem;
  if (MPI_Pack(&offset,              1, MPI_INT,      buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed");
  offset = (ad->csline  == NULL)  ? -1 : ad->csline - ad->mem;
  if (MPI_Pack(&offset,              1, MPI_INT,      buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed");
  offset = (ad->model   == NULL)  ? -1 : ad->model - ad->mem;
  if (MPI_Pack(&offset,              1, MPI_INT,      buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed");
  offset = (ad->mline   == NULL)  ? -1 : ad->mline - ad->mem;
  if (MPI_Pack(&offset,              1, MPI_INT,      buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed");
  offset = (ad->aseq    == NULL)     ? -1 : ad->aseq - ad->mem;
  if (MPI_Pack(&offset,              1, MPI_INT,      buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed");
  offset = (ad->ppline  == NULL)  ? -1 : ad->ppline - ad->mem;
  if (MPI_Pack(&offset,              1, MPI_INT,      buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&ad->N,               1, MPI_INT,      buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed");
  offset = (ad->cmname == NULL)   ? -1 : ad->cmname - ad->mem;
  if (MPI_Pack(&offset,              1, MPI_INT,      buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed");
  offset = (ad->cmacc  == NULL)   ? -1 : ad->cmacc - ad->mem;
  if (MPI_Pack(&offset,              1, MPI_INT,      buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed");
  offset = (ad->cmdesc == NULL)   ? -1 : ad->cmdesc - ad->mem;
  if (MPI_Pack(&offset,              1, MPI_INT,      buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed");

  if (MPI_Pack(&ad->cfrom,           1, MPI_INT,      buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&ad->cto,             1, MPI_INT,      buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&ad->clen,            1, MPI_INT,      buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed");
  offset = (ad->sqname  == NULL)  ? -1 : ad->sqname - ad->mem;
  if (MPI_Pack(&offset,              1, MPI_INT,      buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed");
  offset = (ad->sqacc   == NULL)  ? -1 : ad->sqacc - ad->mem;
  if (MPI_Pack(&offset,              1, MPI_INT,      buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed");
  offset = (ad->sqdesc  == NULL)  ? -1 : ad->sqdesc - ad->mem;
  if (MPI_Pack(&offset,              1, MPI_INT,      buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&ad->sqfrom,          1, MPI_LONG,     buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&ad->sqto,            1, MPI_LONG,     buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&ad->L,               1, MPI_LONG,     buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed");

  if (MPI_Pack(&ad->used_optacc,     1, MPI_INT,      buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&ad->aln_sc,          1, MPI_FLOAT,    buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&ad->used_hbands,     1, MPI_INT,      buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&ad->matrix_Mb,       1, MPI_FLOAT,    buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&ad->elapsed_secs,    1, MPI_DOUBLE,   buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed");

  if (MPI_Pack(&ad->memsize,         1, MPI_INT,      buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack( ad->mem,   ad->memsize, MPI_CHAR,     buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  
  if (*pos > n) ESL_EXCEPTION(eslEMEM, "buffer overflow");
  return eslOK;
  
 ERROR:
  return status;
}

/* Function:  cm_hit_MPIUnpack()
 * Synopsis:  Unpacks an HIT from an MPI buffer.
 * Incept:    EPN, Thu Jun  2 14:29:20 2011 
 *            MSF, Mon Sep  21 22:22:00 2009 [Janelia] (p7_hit_MPIUnpack())
 *
 * Purpose:   Unpack a newly allocated HIT from MPI packed buffer
 *            <buf>, starting from position <*pos>, where the total length
 *            of the buffer in bytes is <n>. 
 *            
 * Returns:   <eslOK> on success. <*pos> is updated to the position of
 *            the next element in <buf> to unpack (if any). <*ret_hit>
 *            contains a newly allocated HIT, which the caller is
 *            responsible for free'ing.
 *            
 * Throws:    <eslESYS> on an MPI call failure. <eslEMEM> on allocation failure.
 *            In either case, <*ret_hit> is <NULL>, and the state of <buf>
 *            and <*pos> is undefined and should be considered to be corrupted.
 */
int
cm_hit_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, CM_HIT *hit)
{
  int  status;
  int  rfline, nline, csline, model, mline, aseq, ppline;
  int  cmname, cmacc, cmdesc;
  int  sqname, sqacc, sqdesc;

  CM_ALIDISPLAY *ad; 
  ESL_ALLOC(ad, sizeof(CM_ALIDISPLAY));

  if (MPI_Unpack(buf, n, pos, &hit->start,       1, MPI_LONG,   comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
  if (MPI_Unpack(buf, n, pos, &hit->stop,        1, MPI_LONG,   comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
  if (MPI_Unpack(buf, n, pos, &hit->in_rc,       1, MPI_INT,    comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
  if (MPI_Unpack(buf, n, pos, &hit->seq_idx,     1, MPI_LONG,   comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
  if (MPI_Unpack(buf, n, pos, &hit->score,       1, MPI_FLOAT,  comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
  if (MPI_Unpack(buf, n, pos, &hit->pvalue,      1, MPI_DOUBLE, comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
  if (MPI_Unpack(buf, n, pos, &hit->evalue,      1, MPI_DOUBLE, comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
  if (MPI_Unpack(buf, n, pos, &hit->flags,       1, MPI_INT,    comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  if ((status = esl_mpi_UnpackOpt(buf, n, pos,   (void**)&(hit->name),        NULL, MPI_CHAR,  comm)) != eslOK) goto ERROR;
  if ((status = esl_mpi_UnpackOpt(buf, n, pos,   (void**)&(hit->acc),         NULL, MPI_CHAR,  comm)) != eslOK) goto ERROR;
  if ((status = esl_mpi_UnpackOpt(buf, n, pos,   (void**)&(hit->desc),        NULL, MPI_CHAR,  comm)) != eslOK) goto ERROR;

  if (MPI_Unpack(buf, n, pos, &rfline,             1, MPI_INT,    comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &nline,              1, MPI_INT,    comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &csline,             1, MPI_INT,    comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &model,              1, MPI_INT,    comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &mline,              1, MPI_INT,    comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &aseq,               1, MPI_INT,    comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &ppline,             1, MPI_INT,    comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &ad->N,              1, MPI_INT,    comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &cmname,             1, MPI_INT,    comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &cmacc,              1, MPI_INT,    comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &cmdesc,             1, MPI_INT,    comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &ad->cfrom,          1, MPI_INT,    comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &ad->cto,            1, MPI_INT,    comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &ad->clen,           1, MPI_INT,    comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &sqname,             1, MPI_INT,    comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &sqacc,              1, MPI_INT,    comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &sqdesc,             1, MPI_INT,    comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &ad->sqfrom,         1, MPI_LONG,   comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &ad->sqto,           1, MPI_LONG,   comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &ad->L,              1, MPI_LONG,   comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  if (MPI_Unpack(buf, n, pos, &ad->used_optacc,    1, MPI_INT,    comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &ad->aln_sc,         1, MPI_FLOAT,  comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &ad->used_hbands,    1, MPI_INT,    comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &ad->matrix_Mb,      1, MPI_FLOAT,  comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &ad->elapsed_secs,   1, MPI_DOUBLE, comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  if (MPI_Unpack(buf, n, pos, &ad->memsize,        1, MPI_INT,    comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  /* allocate the string pools for the alignments */
  ESL_ALLOC(ad->mem, ad->memsize);
  if (MPI_Unpack(buf, n, pos,  ad->mem,  ad->memsize, MPI_CHAR,   comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 

  ad->rfline  = (rfline == -1)  ? NULL : ad->mem + rfline;
  ad->nline   = (nline == -1)   ? NULL : ad->mem + nline;
  ad->csline  = (csline == -1)  ? NULL : ad->mem + csline;
  ad->model   = (model == -1)   ? NULL : ad->mem + model;
  ad->mline   = (mline == -1)   ? NULL : ad->mem + mline;
  ad->aseq    = (aseq == -1)    ? NULL : ad->mem + aseq;
  ad->ppline  = (ppline == -1)  ? NULL : ad->mem + ppline;

  ad->cmname  = (cmname == -1)  ? NULL : ad->mem + cmname;
  ad->cmacc   = (cmacc == -1)   ? NULL : ad->mem + cmacc;
  ad->cmdesc  = (cmdesc == -1)  ? NULL : ad->mem + cmdesc; 

  ad->sqname  = (sqname == -1)  ? NULL : ad->mem + sqname;
  ad->sqacc   = (sqacc == -1)   ? NULL : ad->mem + sqacc;
  ad->sqdesc  = (sqdesc == -1)  ? NULL : ad->mem + sqdesc;

  hit->ad = ad;

  return eslOK;

 ERROR:
  return status;
}

#endif
