/* Optional support for MPI parallelization in Infernal.
 * 
 * Contents:
 *    1. Communicating a CM.
 *    2. Benchmark driver.
 *    3. Unit tests.
 *    4. Test driver.
 *    5. Copyright and license information.
 *    
 * EPN, Mon Aug 27 12:38:13 2007
 * SVN $Id$
 */

#include "esl_config.h"
#include "config.h"

#ifdef HAVE_MPI

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "mpi.h"

#include "easel.h"
#include "esl_mpi.h"

#include "structs.h"
#include "funcs.h"
#include "cm_dispatch.h"
#include "stats.h"

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
    if ((status = cm_MPIPackSize(cm, comm, &sz)) != eslOK) return status;
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
    if ((status = cm_MPIPack(cm, *buf, n, &pos, comm)) != eslOK) return status;
    }

  /* Broadcast the size of the packed CM */
  if (MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD) != 0) ESL_EXCEPTION(eslESYS, "mpi broadcast failed");
  /* Broadcast the packed CM */
  if (MPI_Bcast (*buf, n, MPI_PACKED, 0, MPI_COMM_WORLD) != 0) ESL_EXCEPTION(eslESYS, "mpi broadcast failed");
  return eslOK;
  
 ERROR:
  return status;
}

/* Function:  cm_worker_MPIUnpack()
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

  return cm_MPIUnpack(abc, *buf, *nalloc, &pos, comm, ret_cm);

 ERROR:
  *ret_cm = NULL;
  return status;
}

/* Function:  cm_worker_MPIUnpack()
 * Synopsis:  Unpacks a CM <cm> from an MPI buffer.
 * Incept:    EPN, Mon Aug 27 14:44:11 2007
 *
 * Purpose:   Unpack a newly allocated CM from MPI packed buffer
 *            <buf>, starting from position <*pos>, where the total length
 *            of the buffer in bytes is <n>. 
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
cm_MPIUnpack(ESL_ALPHABET **abc, char *buf, int n, int *pos, MPI_Comm comm, CM_t **ret_cm)
{
  int     status;
  CM_t *cm = NULL;
  int M, nnodes, K, atype;

  if (MPI_Unpack(buf, n, pos, &M,      1, MPI_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &nnodes, 1, MPI_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &atype,  1, MPI_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  /* Set or verify the alphabet */
  if (*abc == NULL)	{	/* still unknown: set it, pass control of it back to caller */
    if ((*abc = esl_alphabet_Create(atype)) == NULL)       { status = eslEMEM;      goto ERROR; }
  } else {			/* already known: check it */
    if ((*abc)->type != atype)                             { status = eslEINCOMPAT; goto ERROR; }
  }

  if ((cm = CreateCM(nnodes,M,(*abc))) == NULL) { status = eslEMEM; goto ERROR;    }
  K = cm->abc->K;

  /* Unpack the rest of the CM */
  if (MPI_Unpack(buf, n, pos, &(cm->flags),              1,   MPI_INT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->config_opts),        1,   MPI_INT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->search_opts),        1,   MPI_INT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->align_opts),         1,   MPI_INT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->nseq),               1,   MPI_INT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->clen),               1,   MPI_INT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->iel_selfsc),         1,   MPI_INT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->W),                  1,   MPI_INT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->enf_start),          1,   MPI_INT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->cutoff_type),        1,   MPI_INT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->cp9_cutoff_type),    1,   MPI_INT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->hmmpad),             1,   MPI_INT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");

  if (MPI_Unpack(buf, n, pos, &(cm->el_selfsc),          1, MPI_FLOAT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->enf_scdiff),         1, MPI_FLOAT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->sc_boost),           1, MPI_FLOAT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->cp9_sc_boost),       1, MPI_FLOAT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->cutoff),             1, MPI_FLOAT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->cp9_cutoff),         1, MPI_FLOAT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->pbegin),             1, MPI_FLOAT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->pend),               1, MPI_FLOAT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->eff_nseq),           1, MPI_FLOAT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->ga),                 1, MPI_FLOAT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->tc),                 1, MPI_FLOAT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->nc),                 1, MPI_FLOAT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, cm->null,                  K, MPI_FLOAT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");

  if (MPI_Unpack(buf, n, pos, &(cm->beta),               1,MPI_DOUBLE, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->tau),                1,MPI_DOUBLE, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");

  if (MPI_Unpack(buf, n, pos, cm->sttype,            (M+1),  MPI_CHAR, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, cm->stid,              (M+1),  MPI_CHAR, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");

  if (MPI_Unpack(buf, n, pos, cm->ndidx,                 M,   MPI_INT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, cm->cfirst,                M,   MPI_INT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, cm->cnum,                  M,   MPI_INT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, cm->plast,                 M,   MPI_INT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, cm->pnum,                  M,   MPI_INT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, cm->ibeginsc,              M,   MPI_INT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, cm->iendsc,                M,   MPI_INT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");

  if (MPI_Unpack(buf, n, pos, cm->begin,                 M, MPI_FLOAT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, cm->end,                   M, MPI_FLOAT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, cm->beginsc,               M, MPI_FLOAT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, cm->endsc,                 M, MPI_FLOAT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  
  if (MPI_Unpack(buf, n, pos, cm->nodemap,          nnodes,   MPI_INT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, cm->ndtype,           nnodes,   MPI_INT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  
  if (MPI_Unpack(buf, n, pos, cm->e[0],              M*K*K, MPI_FLOAT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, cm->esc[0],            M*K*K, MPI_FLOAT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");

  if (MPI_Unpack(buf, n, pos, cm->t[0],       M*MAXCONNECT, MPI_FLOAT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, cm->tsc[0],     M*MAXCONNECT, MPI_FLOAT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  
  if (MPI_Unpack(buf, n, pos, cm->iesc[0],           M*K*K,   MPI_INT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");

  if (MPI_Unpack(buf, n, pos, cm->itsc[0],    M*MAXCONNECT,   MPI_INT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");

  if (                            (status = esl_mpi_UnpackOpt(buf, n, pos, (void**)&(cm->enf_seq), NULL, MPI_CHAR, comm)) != eslOK) goto ERROR;
  if (                            (status = esl_mpi_UnpackOpt(buf, n, pos, (void**)&(cm->name),    NULL, MPI_CHAR, comm)) != eslOK) goto ERROR;
  if (cm->flags & CMH_ACC)  { if ((status = esl_mpi_UnpackOpt(buf, n, pos, (void**)&(cm->acc),     NULL, MPI_CHAR, comm)) != eslOK) goto ERROR; }
  if (cm->flags & CMH_DESC) { if ((status = esl_mpi_UnpackOpt(buf, n, pos, (void**)&(cm->desc),    NULL, MPI_CHAR, comm)) != eslOK) goto ERROR; }
  
  
  *ret_cm = cm;
  return eslOK;
  
  ERROR:
  if (cm != NULL) FreeCM(cm);
  return status;
}


/* Function:  cm_MPIPack()
 * Incept:    EPN, Mon Aug 27 13:57:54 2007
 *
 * Purpose:   Packs CM <cm> into an MPI packed message buffer <buf>
 *            of length <n> bytes, starting at byte position <*position>,
 *            for MPI communicator <comm>.
 *            
 *            The caller must know that <buf>'s allocation of <n>
 *            bytes is large enough to append the packed CM at
 *            position <*pos>. This typically requires a call to
 *            <cm_MPIPackSize()> first, and reallocation if
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
cm_MPIPack(CM_t *cm, char *buf, int n, int *pos, MPI_Comm comm)
{
  int   status;
  int   K      = cm->abc->K;
  int   M      = cm->M;
  int   nnodes = cm->nodes;

  if (MPI_Pack(&M,                        1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&nnodes,                   1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack((void *) &(cm->abc->type), 1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->flags),              1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->config_opts),        1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->search_opts),        1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->align_opts),         1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->nseq),               1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->clen),               1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->iel_selfsc),         1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->W),                  1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->enf_start),          1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->cutoff_type),        1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->cp9_cutoff_type),    1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->hmmpad),             1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");

  if (MPI_Pack(&(cm->el_selfsc),          1, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->enf_scdiff),         1, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->sc_boost),           1, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->cp9_sc_boost),       1, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->cutoff),             1, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->cp9_cutoff),         1, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->pbegin),             1, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->pend),               1, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->eff_nseq),           1, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->ga),                 1, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->tc),                 1, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->nc),                 1, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(cm->null,                  K, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");

  if (MPI_Pack(&(cm->beta),               1,MPI_DOUBLE, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->tau),                1,MPI_DOUBLE, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");

  if (MPI_Pack(cm->sttype,            (M+1),  MPI_CHAR, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(cm->stid,              (M+1),  MPI_CHAR, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");

  if (MPI_Pack(cm->ndidx,                 M,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(cm->cfirst,                M,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(cm->cnum,                  M,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(cm->plast,                 M,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(cm->pnum,                  M,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(cm->ibeginsc,              M,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(cm->iendsc,                M,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");

  if (MPI_Pack(cm->begin,                 M, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(cm->end,                   M, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(cm->beginsc,               M, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(cm->endsc,                 M, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");

  if (MPI_Pack(cm->nodemap,          nnodes,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(cm->ndtype,           nnodes,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");

  if (MPI_Pack(cm->e[0],              M*K*K, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(cm->esc[0],            M*K*K, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");

  if (MPI_Pack(cm->t[0],       M*MAXCONNECT, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(cm->tsc[0],     M*MAXCONNECT, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  
  if (MPI_Pack(cm->iesc[0],           M*K*K,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");

  if (MPI_Pack(cm->itsc[0],    M*MAXCONNECT,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");

  
  if (                            (status = esl_mpi_PackOpt(cm->enf_seq, -1, MPI_CHAR, buf, n, pos, comm)) != eslOK) return status; 
  if (                            (status = esl_mpi_PackOpt(cm->name,    -1, MPI_CHAR, buf, n, pos, comm)) != eslOK) return status; 
  if (cm->flags & CMH_ACC)  { if ((status = esl_mpi_PackOpt(cm->acc,     -1, MPI_CHAR, buf, n, pos, comm)) != eslOK) return status; }
  if (cm->flags & CMH_DESC) { if ((status = esl_mpi_PackOpt(cm->desc,    -1, MPI_CHAR, buf, n, pos, comm)) != eslOK) return status; }

  if (*pos > n) ESL_EXCEPTION(eslEMEM, "buffer overflow");
  return eslOK;
}

/* Function:  cm_justread_MPIPack()
 * Incept:    EPN, Mon Aug 27 14:24:57 2007
 *
 * Purpose:   Packs CM <cm> into an MPI packed message buffer <buf>
 *            of length <n> bytes, starting at byte position <*position>,
 *            for MPI communicator <comm>. Differs from cm_MPIPack()
 *            in that the full CM data structure is not packed, only
 *            the parts that could have been changed from their initial
 *            (default) values in CMFileRead() are packed.
 *            
 *            The caller must know that <buf>'s allocation of <n>
 *            bytes is large enough to append the packed CM at
 *            position <*pos>. This typically requires a call to
 *            <cm_MPIPackSize()> first, and reallocation if
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
  int   K      = cm->abc->K;
  int   M      = cm->M;
  int   nnodes = cm->nodes;

  if (MPI_Pack(&M,                        1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&nnodes,                   1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack((void *) &(cm->abc->type), 1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->nseq),               1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->clen),               1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");

  if (MPI_Pack(&(cm->el_selfsc),          1, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->eff_nseq),           1, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->ga),                 1, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->tc),                 1, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->nc),                 1, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(cm->null,                  K, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");

  if (MPI_Pack(cm->sttype,            (M+1),  MPI_CHAR, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(cm->stid,              (M+1),  MPI_CHAR, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");

  if (MPI_Pack(cm->ndidx,                 M,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(cm->cfirst,                M,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(cm->cnum,                  M,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(cm->plast,                 M,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(cm->pnum,                  M,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");

  if (MPI_Pack(cm->nodemap,          nnodes,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(cm->ndtype,           nnodes,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");

  if (MPI_Pack(cm->e[0],              M*K*K, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");

  if (MPI_Pack(cm->t[0],       M*MAXCONNECT, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  
  if (                            (status = esl_mpi_PackOpt(cm->name,    -1, MPI_CHAR, buf, n, pos, comm)) != eslOK) return status; 
  if (cm->flags & CMH_ACC)  { if ((status = esl_mpi_PackOpt(cm->acc,     -1, MPI_CHAR, buf, n, pos, comm)) != eslOK) return status; }
  if (cm->flags & CMH_DESC) { if ((status = esl_mpi_PackOpt(cm->desc,    -1, MPI_CHAR, buf, n, pos, comm)) != eslOK) return status; }

  if (*pos > n) ESL_EXCEPTION(eslEMEM, "buffer overflow");
  return eslOK;
}

/* Function:  cm_MPIPackSize()
 * Synopsis:  Calculates size needed to pack a CM.
 * Incept:    EPN, Mon Aug 27 10:34:15 2007
 *            based on p7_hmm_MPIPackSize() from HMMER3.
 *
 * Purpose:   Calculate an upper bound on the number of bytes
 *            that <cm_MPIPack()> will need to pack a CM
 *            <cm> in a packed MPI message for MPI communicator
 *            <comm>; return that number of bytes in <*ret_n>.
 *
 * Returns:   <eslOK> on success, and <*ret_n> contains the answer.
 *
 * Throws:    <eslESYS> if an MPI call fails, and <*ret_n> is 0.
 */
int
cm_MPIPackSize(CM_t *cm, MPI_Comm comm, int *ret_n)
{
  int   status;
  int   n = 0;
  int   K = cm->abc->K;
  int   M = cm->M;
  int   nnodes = cm->nodes;
  int   sz;

  if (MPI_Pack_size(1,         MPI_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");   
  n += 15*sz; 
  /* M, nodes, abc->type, flags, config_opts, search_opts, align_opts, nseq, clen, iel_selfsc, (10)
   * W, enf_start, cutoff_type, cp9_cutoff_type, hmmpad, (5) 
   */

  if (MPI_Pack_size(1,       MPI_FLOAT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += (12+K)*sz; 
  /* el_selfsc, enf_scdiff, sc_boost, cp9_sc_boost, cutoff, cp9_cutoff, pbegin, pend, eff_nseq (9) 
   * ga, tc, nc, null (3+K) */

  if (MPI_Pack_size(1,      MPI_DOUBLE, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += 2*sz; 
  /* beta, tau */

  if (MPI_Pack_size((M+1),        MPI_CHAR, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += 2*sz; 
  /* sttype, stid */

  if (MPI_Pack_size(M,         MPI_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += 7*sz; 
  /* ndidx, cfirst, cnum, plast, pnum, ibeginsc, iendsc */

  if (MPI_Pack_size(M,       MPI_FLOAT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += 4*sz; 
  /* begin, end, beginsc, endsc */

  if (MPI_Pack_size(nnodes,    MPI_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += 2*sz; 
  /* nodemap, ndtype */

  if (MPI_Pack_size(M*K*K,   MPI_FLOAT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += 2*sz; 
  /* e, esc */

  if (MPI_Pack_size(M*MAXCONNECT,MPI_FLOAT,comm,&sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += 2*sz; 
  /* t, tsc */

  if (MPI_Pack_size(M*K*K,    MPI_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += sz; 
  /* iesc */

  if (MPI_Pack_size(M*MAXCONNECT,MPI_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += sz; 
  /* itsc */

  if ((status = esl_mpi_PackOptSize(cm->enf_seq, -1, MPI_CHAR, comm, &sz)) != eslOK) goto ERROR; 
  n += sz; /* enf_seq */

  if ((status = esl_mpi_PackOptSize(cm->name, -1, MPI_CHAR, comm, &sz)) != eslOK) goto ERROR; 
  n += sz; /* name */

  if (cm->flags & CMH_ACC)  { if ((status = esl_mpi_PackOptSize(cm->acc, -1,  MPI_CHAR,  comm, &sz)) != eslOK) goto ERROR;  n+= sz; }
  if (cm->flags & CMH_DESC) { if ((status = esl_mpi_PackOptSize(cm->desc,-1,  MPI_CHAR,  comm, &sz)) != eslOK) goto ERROR;  n+= sz; }
  
  *ret_n = n;
  return eslOK;

 ERROR:
  *ret_n = 0;
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
 *            that <cm_MPIPack()> will need to pack a CM
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
  int   n = 0;
  int   K = cm->abc->K;
  int   M = cm->M;
  int   nnodes = cm->nodes;
  int   sz;

  if (MPI_Pack_size(1,         MPI_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");   
  n += 5*sz; /* M, nodes, abc->type, nseq, clen */ 

  if (MPI_Pack_size(1,       MPI_FLOAT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += (5+K)*sz; 
  /* el_selfsc, eff_nseq, ga, tc, nc, null (5+K) */

  if (MPI_Pack_size(M,        MPI_CHAR, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += 2*sz; 
  /* sttype, stid */

  if (MPI_Pack_size(M,         MPI_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += 5*sz; 
  /* ndidx, cfirst, cnum, plast, pnum, */

  if (MPI_Pack_size(nnodes,    MPI_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += 2*sz; 
  /* nodemap, ndtype */

  if (MPI_Pack_size(M*K*K,   MPI_FLOAT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += sz; 
  /* e */

  if (MPI_Pack_size(M*MAXCONNECT,MPI_FLOAT,comm,&sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += sz; 
  /* t */

  if ((status = esl_mpi_PackOptSize(cm->name, -1, MPI_CHAR, comm, &sz)) != eslOK) goto ERROR; 
  n += sz; /* name */

  if (cm->flags & CMH_ACC)  { if ((status = esl_mpi_PackOptSize(cm->acc, -1,  MPI_CHAR,  comm, &sz)) != eslOK) goto ERROR;  n+= sz; }
  if (cm->flags & CMH_DESC) { if ((status = esl_mpi_PackOptSize(cm->desc,-1,  MPI_CHAR,  comm, &sz)) != eslOK) goto ERROR;  n+= sz; }
  
  *ret_n = n;
  return eslOK;

 ERROR:
  *ret_n = 0;
  return status;

}

/**************************************************************************************/
/* EPN, Thu May 10 10:11:18 2007 New functions roughly following Easel/H3 conventions */
/* Function: mpi_worker_search_target()
 * Incept:   EPN, Wed May  9 17:07:48 2007
 * Purpose:  The main control for an MPI worker process for searching sequences. 
 *           Worker receives CM, then loops over receipt of sequences, returning
 *           best score and results data structure for each.
 *           Never do revcomp, we'll call this function twice once with 
 *           plus once with minus strand.
 */
void
mpi_worker_search_target(CM_t *cm, int my_rank)
{
  int status;
  char *dsq = NULL;
  int   L;
  float best_sc;

  int doing_cm_stats  = FALSE;
  int doing_cp9_stats = FALSE;

  int ctr = 0;
  if(cm->search_opts & CM_SEARCH_HMMONLY) doing_cp9_stats = TRUE;
  else doing_cm_stats = TRUE;
  /* Main loop */
  while (dsq_MPIRecv(&dsq, &L) == eslOK)
    {
      best_sc = actually_search_target(cm, dsq, 1, L, 
				       0.,    /* minimum CM bit cutoff, irrelevant (?) */
				       0.,    /* minimum CP9 bit cutoff, irrelevant (?) */
				       NULL,  /* do not keep results */
				       FALSE, /* do not filter with a CP9 HMM */
				       doing_cm_stats, doing_cp9_stats,
				       NULL); /* filter fraction, nobody cares */
      
      MPI_Send(&(best_sc), 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
      free(dsq);
    }
  return;

 ERROR:
  if (dsq != NULL) free(dsq);
  return;
}

/* Function:  dsq_MPISend()
 * Incept:    EPN, Wed May  9 17:30:14 2007
 *
 * Purpose:   Send sequence <dsq> to processor <dest>.
 *            
 *            If <dsq> is NULL, sends a end-of-data signal to <dest>, to
 *            tell it to shut down.
 */
int
dsq_MPISend(char *dsq, int L, int dest)
{
  if(dsq == NULL) 
    {
      int eod = -1;
      MPI_Send(&eod, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
      return eslOK;
    }
  MPI_Send(&(L), 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
  /* receiver will now allocate storage, before reading on...*/
  MPI_Send(dsq, (L+2), MPI_CHAR, dest, 0, MPI_COMM_WORLD);
  return eslOK;
}

/* Function:  dsq_MPISend()
 * Incept:    EPN, Thu Jun  7 15:02:34 2007
 *
 * Purpose:   Send sequence <dsq> and max bit sc <maxsc> to processor <dest>.
 *            
 *            If <dsq> is NULL, sends a end-of-data signal to <dest>, to
 *            tell it to shut down.
 */
int
dsq_maxsc_MPISend(char *dsq, int L, float maxsc, int dest)
{
  if(dsq == NULL) 
    {
      int eod = -1;
      MPI_Send(&eod, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
      return eslOK;
    }
  MPI_Send(&(L), 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
  /* receiver will now allocate storage, before reading on...*/
  MPI_Send(dsq, (L+2), MPI_CHAR, dest, 0, MPI_COMM_WORLD);
  MPI_Send(&(maxsc), 1, MPI_FLOAT, dest, 0, MPI_COMM_WORLD);
  return eslOK;
}

/* Function:  dsq_MPIRecv()
 * Incept:    EPN, Wed May  9 17:34:43 2007
 *
 * Purpose:   Receive a sequence sent from the master MPI process (src=0)
 *            on a worker MPI process. 
 *            
 *            If it receives an end-of-data signal, returns <eslEOD>.
 */
int
dsq_MPIRecv(char **ret_dsq, int *ret_L)
{
  int status;
  char *dsq = NULL;
  MPI_Status mpistatus;
  int L;
  
  MPI_Recv(&L, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpistatus);
  if (L == -1) return eslEOD;
  ESL_ALLOC(dsq, sizeof(char) * (L+2));
  MPI_Recv(dsq, (L+2), MPI_CHAR, 0, 0, MPI_COMM_WORLD, &mpistatus);
  *ret_L   = L;
  *ret_dsq = dsq;
  return eslOK;

 ERROR:
  return status;
}

/* Function:  dsq_maxsc_MPIRecv()
 * Incept:    EPN, Thu Jun  7 15:00:29 2007    
 *
 * Purpose:   Receive a sequence and maximum score 
 *            sent from the master MPI process (src=0)
 *            on a worker MPI process. 
 *            
 *            If it receives an end-of-data signal, returns <eslEOD>.
 */
int
dsq_maxsc_MPIRecv(char **ret_dsq, int *ret_L, float *ret_maxsc)
{
  int status;
  char *dsq = NULL;
  MPI_Status mpistatus;
  int L;
  float maxsc;
  
  MPI_Recv(&L, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpistatus);
  if (L == -1) return eslEOD;
  ESL_ALLOC(dsq, sizeof(char) * (L+2));
  MPI_Recv(dsq, (L+2), MPI_CHAR, 0, 0, MPI_COMM_WORLD, &mpistatus);
  MPI_Recv(&maxsc, 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &mpistatus);
  *ret_L   = L;
  *ret_dsq = dsq;
  *ret_maxsc = maxsc;
  return eslOK;

 ERROR:
  return status;
}

/* Function: mpi_worker_cm_and_cp9_search()
 * Incept:   EPN, Thu May 10 10:04:02 2007
 * Purpose:  The main control for an MPI worker process for searching sequences
 *           twice, once with a CM and once with a CP9, both scores are returned.
 *           Called in mpi_FindCP9FilterThreshold9).
 * Args:
 *           cm       - the covariance model
 *           do_fast  - don't search with CM, only do CP9 search
 *           my_rank  - my MPI rank
 */
void
mpi_worker_cm_and_cp9_search(CM_t *cm, int do_fast, int my_rank)
{
  int status;
  char *dsq = NULL;
  int   L;
  float *scores = NULL;
  ESL_ALLOC(scores, sizeof(float) * 2);
  int was_hmmonly;
  if(cm->search_opts & CM_SEARCH_HMMONLY) was_hmmonly = TRUE;
  else was_hmmonly = FALSE;
  /* Main loop */
  while (dsq_MPIRecv(&dsq, &L) == eslOK)
    {
      /* Do the CM search first */
      cm->search_opts &= ~CM_SEARCH_HMMONLY;
      if(do_fast)
	scores[0] = 0.;
      else
	scores[0] = actually_search_target(cm, dsq, 1, L, 
					   0.,    /* minimum CM bit cutoff, irrelevant (?) */
					   0.,    /* minimum CP9 bit cutoff, irrelevant (?) */
					   NULL,  /* do not keep results */
					   FALSE, /* do not filter with a CP9 HMM */
					   FALSE, FALSE, /* not doing CM or CP9 Gumbel calcs */
					   NULL); /* filter fraction, nobody cares */
      /* DO NOT CALL actually_search_target b/c that will run Forward then 
       * Backward to get score of best hit, but we'll be detecting by a
       * Forward scan (then running Backward only on hits above our threshold),
       * since we're calc'ing the threshold here it's impt we only do Forward.
       */
      scores[1] =  CP9Forward(cm, dsq, 1, L, cm->W, 0., 
			      NULL,   /* don't return scores of hits */
			      NULL,   /* don't return posns of hits */
			      NULL,   /* don't keep track of hits */
			      TRUE,   /* we're scanning */
			      FALSE,  /* we're not ultimately aligning */
			      FALSE,  /* we're not rescanning */
			      TRUE,   /* be memory efficient */
			      NULL);  /* don't want the DP matrix back */
      
      MPI_Send(scores, 2, MPI_FLOAT, 0, 0, MPI_COMM_WORLD); /* send together so results don't interleave */
      free(dsq);
    }
  if(was_hmmonly) cm->search_opts |= CM_SEARCH_HMMONLY;
  else cm->search_opts &= ~CM_SEARCH_HMMONLY;
  if(scores != NULL) free(scores);
  /*printf("\trank: %d RETURNING!\n", my_rank);*/
  return;

 ERROR:
  if (dsq != NULL) free(dsq);
  if (scores != NULL) free(scores);
  return;
}

/* Function: mpi_worker_cm_and_cp9_search_maxsc()
 * Incept:   EPN, Thu Jun  7 15:02:54 2007   
 * Purpose:  The main control for an MPI worker process for searching sequences
 *           with decreasingly fast techniques, quitting if any technique 
 *           returns a score greater than some specified maximum bit score. 
 *           The goal is to determine if the optimal parse is within a 
 *           given range during a empirical HMM filter threshold calculation. 
 *           Called in mpi_FindCP9FilterThreshold().
 * Args:
 *           cm       - the covariance model
 *           do_fast  - don't search with CM, only do CP9 search
 *           my_rank  - my MPI rank
 */
void
mpi_worker_cm_and_cp9_search_maxsc(CM_t *cm, int do_fast, int do_minmax, int my_rank)
{
  int status;
  char *dsq = NULL;
  int   L;
  float maxsc;
  float *scores = NULL;
  ESL_ALLOC(scores, sizeof(float) * 2);
  int was_hmmonly;
  int was_hbanded;
  float orig_tau;
  int passed_flag = FALSE;
  orig_tau = cm->tau;

  if(cm->search_opts & CM_SEARCH_HMMONLY) was_hmmonly = TRUE;
  else was_hmmonly = FALSE;
  if(cm->search_opts & CM_SEARCH_HBANDED) was_hbanded = TRUE;
  else was_hbanded = FALSE;
  /* Main loop */
  while (dsq_maxsc_MPIRecv(&dsq, &L, &maxsc) == eslOK)
    {
      /* Do the CM search first */
      cm->search_opts &= ~CM_SEARCH_HMMONLY;
      if(do_fast)
	scores[0] = 0.;
      else if(do_minmax)
	{
	  cm->search_opts |= CM_SEARCH_HBANDED;
	  cm->tau = 0.1;
	  scores[0] = actually_search_target(cm, dsq, 1, L,
					 0.,    /* cutoff is 0 bits (actually we'll find highest
						 * negative score if it's < 0.0) */
					 0.,    /* CP9 cutoff is 0 bits */
					 NULL,  /* don't keep results */
					 FALSE, /* don't filter with a CP9 HMM */
					 FALSE, /* we're not calcing CM  stats */
					 FALSE, /* we're not calcing CP9 stats */
					 NULL); /* filter fraction N/A */
	  if(scores[0] < maxsc) /* search with another, less strict (lower tau)  HMM banded parse */
	    {
	      cm->tau = 1e-10;
	      scores[0] = actually_search_target(cm, dsq, 1, L,
					     0.,    /* cutoff is 0 bits (actually we'll find highest
						     * negative score if it's < 0.0) */
					     0.,    /* CP9 cutoff is 0 bits */
					     NULL,  /* don't keep results */
					     FALSE, /* don't filter with a CP9 HMM */
					     FALSE, /* we're not calcing CM  stats */
					     FALSE, /* we're not calcing CP9 stats */
					     NULL); /* filter fraction N/A */
	    }
	}
      else
	scores[0] = actually_search_target(cm, dsq, 1, L, 
					   0.,    /* minimum CM bit cutoff, irrelevant (?) */
					   0.,    /* minimum CP9 bit cutoff, irrelevant (?) */
					   NULL,  /* do not keep results */
					   FALSE, /* do not filter with a CP9 HMM */
					   FALSE, FALSE, /* not doing CM or CP9 Gumbel calcs */
					   NULL); /* filter fraction, nobody cares */

      /* Now do HMM search, but if do_minmax, only do HMM search 
       * if our CM score hasn't exceeded the max */
      if(do_minmax && scores[0] > maxsc)
	scores[1] = 0.;
      else
	/* DO NOT CALL actually_search_target b/c that will run Forward then 
	 * Backward to get score of best hit, but we'll be detecting by a
	 * Forward scan (then running Backward only on hits above our threshold),
	 * since we're calc'ing the threshold here it's impt we only do Forward.
	 */
	scores[1] =  CP9Forward(cm, dsq, 1, L, cm->W, 0., 
				NULL,   /* don't return scores of hits */
				NULL,   /* don't return posns of hits */
				NULL,   /* don't keep track of hits */
				TRUE,   /* we're scanning */
				FALSE,  /* we're not ultimately aligning */
				FALSE,  /* we're not rescanning */
				TRUE,   /* be memory efficient */
				NULL);  /* don't want the DP matrix back */
      MPI_Send(scores, 2, MPI_FLOAT, 0, 0, MPI_COMM_WORLD); /* send together so results don't interleave */
      free(dsq);
    }
  if(was_hmmonly) cm->search_opts |= CM_SEARCH_HMMONLY;
  else cm->search_opts &= ~CM_SEARCH_HMMONLY;
  if(was_hbanded) cm->search_opts |= CM_SEARCH_HBANDED;
  else cm->search_opts &= ~CM_SEARCH_HBANDED;

  if(scores != NULL) free(scores);
  /*printf("\trank: %d RETURNING!\n", my_rank);*/
  return;
  
 ERROR:
  if (dsq != NULL) free(dsq);
  if (scores != NULL) free(scores);
  return;
}
#endif

/******************************************************************************/
#if 0
void broadcast_cm (CM_t **cm, int mpi_my_rank, int mpi_master_rank) 
{
  char buf[BUFSIZE];      /* Buffer for packing it all but the bulk of the CM */
  int position = 0;         /* Where I am in the buffer */
  int nstates, nnodes;
  int enf_len;
  int nparts;
  int i;
  int p;

  position = 0;
  if (mpi_my_rank == mpi_master_rank) 
    {   /* I'm in charge */
      /* contract check, if we claim to have Gumbel stats, we better have them */
      if((*cm)->flags & CM_GUMBEL_STATS && (*cm)->stats == NULL)
	esl_fatal("ERROR in broadcast_cm() master node claims to have Gumbel stats but cm->stats is NULL!\n");
      nstates = (*cm)->M;
      nnodes = (*cm)->nodes;
      
      /* Basics of the model */
      MPI_Pack (&nstates,                  1, MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD); 
      MPI_Pack (&nnodes,                   1, MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->flags),           1, MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->config_opts),     1, MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->align_opts),      1, MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->search_opts),     1, MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->el_selfsc),       1, MPI_FLOAT,  buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->iel_selfsc),      1, MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->W),               1, MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->enf_start),       1, MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->enf_scdiff),      1, MPI_FLOAT,  buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->sc_boost),        1, MPI_FLOAT,  buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->cp9_sc_boost),    1, MPI_FLOAT,  buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->cutoff_type),     1, MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->cutoff),          1, MPI_FLOAT,  buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->cp9_cutoff_type), 1, MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->cp9_cutoff),      1, MPI_FLOAT,  buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->beta),            1, MPI_DOUBLE, buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->tau),             1, MPI_DOUBLE, buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->hmmpad),          1, MPI_INT,    buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->pbegin),          1, MPI_FLOAT,  buf, BUFSIZE, &position, MPI_COMM_WORLD);
      MPI_Pack (&((*cm)->pend),            1, MPI_FLOAT,  buf, BUFSIZE, &position, MPI_COMM_WORLD);

      /* Take special care with enf_len, this is used later to get cm->enf_seq if nec */
      if((*cm)->enf_start != 0) enf_len = strlen((*cm)->enf_seq);
      else enf_len = 0;
      MPI_Pack (&enf_len,                  1, MPI_INT, buf, BUFSIZE, &position, MPI_COMM_WORLD);

      /* Take special care with number of partitions, used later to get cm->stats if nec */
      if((*cm)->flags & CM_GUMBEL_STATS) nparts = (*cm)->stats->np;
      else nparts = 0;
      MPI_Pack (&nparts,                  1, MPI_INT, buf, BUFSIZE, &position, MPI_COMM_WORLD);

    }
  /* Broadcast to everyone */
  MPI_Bcast (buf, BUFSIZE, MPI_PACKED, mpi_master_rank, MPI_COMM_WORLD);

  /* Decode this first set */
  position = 0;
  if (mpi_my_rank != mpi_master_rank) 
    {
      MPI_Unpack (buf, BUFSIZE, &position, &nstates, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &nnodes, 1, MPI_INT, MPI_COMM_WORLD);
      *cm = CreateCM (nnodes, nstates);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->flags),           1, MPI_INT,   MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->config_opts),     1, MPI_INT,   MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->align_opts),      1, MPI_INT,   MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->search_opts),     1, MPI_INT,   MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->el_selfsc),       1, MPI_FLOAT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->iel_selfsc),      1, MPI_INT,   MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->W),               1, MPI_INT,   MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->enf_start),       1, MPI_INT,   MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->enf_scdiff),      1, MPI_FLOAT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->sc_boost),        1, MPI_FLOAT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->cp9_sc_boost),    1, MPI_FLOAT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->ffract),          1, MPI_FLOAT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->cutoff_type),     1, MPI_INT,   MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->cutoff),          1, MPI_FLOAT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->cp9_cutoff_type), 1, MPI_FLOAT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->cp9_cutoff),      1, MPI_FLOAT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->beta),            1, MPI_DOUBLE,MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->tau),             1, MPI_DOUBLE,MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->hmmpad),          1, MPI_INT,   MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->pbegin),          1, MPI_FLOAT, MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &((*cm)->pend),            1, MPI_FLOAT,   MPI_COMM_WORLD);

      MPI_Unpack (buf, BUFSIZE, &position, &enf_len,                  1, MPI_INT,   MPI_COMM_WORLD);
      MPI_Unpack (buf, BUFSIZE, &position, &nparts,                   1, MPI_INT,   MPI_COMM_WORLD);
    }
  /* Now we broadcast the rest of the model using many calls to MPI_Bcast.  
     This is inefficient, but is probably negligible compared to the actual 
     searches */
  MPI_Bcast ((*cm)->null,   Alphabet_size, MPI_FLOAT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->sttype, nstates,       MPI_CHAR,  mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->ndidx,  nstates,       MPI_INT,   mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->stid,   nstates,       MPI_CHAR,  mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->cfirst, nstates,       MPI_INT,   mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->cnum,   nstates,       MPI_INT,   mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->plast,  nstates,       MPI_INT,   mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->pnum,   nstates,       MPI_INT,   mpi_master_rank, MPI_COMM_WORLD);
 
  MPI_Bcast ((*cm)->nodemap, nnodes, MPI_INT,  mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->ndtype,  nnodes, MPI_CHAR, mpi_master_rank, MPI_COMM_WORLD);

  MPI_Bcast ((*cm)->begin,    nstates, MPI_FLOAT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->end,      nstates, MPI_FLOAT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->beginsc,  nstates, MPI_FLOAT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->endsc,    nstates, MPI_FLOAT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->ibeginsc, nstates, MPI_INT,   mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->iendsc,   nstates, MPI_INT,   mpi_master_rank, MPI_COMM_WORLD);

  /* These next calls depend on Sean's FMX2Alloc to be what CreateCM calls, and to allocate one large
     memory chunk at x[0] (where x is float **) and then fill in x[1]..x[n] with the appropriate offsets into
     this chunk of memory */
  MPI_Bcast ((*cm)->t[0], nstates*MAXCONNECT, MPI_FLOAT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->e[0], nstates*Alphabet_size*Alphabet_size, MPI_FLOAT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->tsc[0], nstates*MAXCONNECT, MPI_FLOAT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->esc[0], nstates*Alphabet_size*Alphabet_size, MPI_FLOAT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->itsc[0], nstates*MAXCONNECT, MPI_INT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast ((*cm)->iesc[0], nstates*Alphabet_size*Alphabet_size, MPI_INT, mpi_master_rank, MPI_COMM_WORLD);

  /* Broadcast the enf_seq, if it's NULL (enf_start == 0) we don't */
  if((*cm)->enf_start != 0)
    {
      if (mpi_my_rank != mpi_master_rank) 
	ESL_ALLOC((*cm)->enf_seq, sizeof(char) * (enf_len+1));
      MPI_Bcast((*cm)->enf_seq, enf_len, MPI_CHAR, mpi_master_rank, MPI_COMM_WORLD);
      if (mpi_my_rank != mpi_master_rank) 
	(*cm)->enf_seq[enf_len] = '\0';
    }

  /* Broadcast the Gumbel stats if they exist 
   * IMPT: currently filter threshold stats are NOT broadcasted as they're only
   * used to get cm->cp9_cutoff, which is broadcasted separately. We could get 
   * away with not broadcasting these stats too - though we'd have to modify 
   * parallel_search_database() to be independent on Gumbel params */
  if((*cm)->flags & CM_GUMBEL_STATS) /* flags were already sent/received */
    {
      if (mpi_my_rank != mpi_master_rank) 
	(*cm)->stats = AllocCMStats(nparts); /* nparts was already sent/recd */
      for(i = 0; i < NGUMBELMODES; i++)
	for(p = 0; p < nparts; p++)
	  {	    
	    MPI_Bcast(&((*cm)->stats->gumAA[i][p]->N),      1, MPI_INT,    mpi_master_rank, MPI_COMM_WORLD);
	    MPI_Bcast(&((*cm)->stats->gumAA[i][p]->L),      1, MPI_INT,    mpi_master_rank, MPI_COMM_WORLD);
	    MPI_Bcast(&((*cm)->stats->gumAA[i][p]->mu),     1, MPI_DOUBLE, mpi_master_rank, MPI_COMM_WORLD);
	    MPI_Bcast(&((*cm)->stats->gumAA[i][p]->lambda), 1, MPI_DOUBLE, mpi_master_rank, MPI_COMM_WORLD);
	  }
    }
  return eslOK;
}  

/*
 * Function: search_results_MPIUnpack() 
 * Date:     RJK, Wed May 29, 2002 [St. Louis]
 * Purpose:  Does a blocking call to MPI_Recv until a process finishes, then
 *           processes results and returns.
 */
int search_check_results (db_seq_t **active_seqs, job_t **process_status, int D) {
  char *buf;
  int bufsize;
  MPI_Status status;
  int data_from;        /* Who's sending us data */
  char results_type;
  int position = 0;
  int num_results;
  int start, stop, bestr;
  float score;
  db_seq_t *cur_seq;
  char in_revcomp;
  int index;
  int i;
  int cur_seq_index;
  Parsetree_t *tr;

  /* Get the size of the buffer */
  MPI_Recv (&bufsize, 1, MPI_INT, MPI_ANY_SOURCE, 
	    SEARCH_STD_SCAN_RESULTS_SIZE_TAG, MPI_COMM_WORLD, &status);
  data_from = status.MPI_SOURCE;
  ESL_ALLOC(buf, (char)*bufsize);

  /* Figure out the sequence it belongs to */
  cur_seq_index = process_status[data_from]->db_seq_index;
  cur_seq = active_seqs[cur_seq_index];
  index = process_status[data_from]->index;
  in_revcomp = process_status[data_from]->in_revcomp;

  /* Clear this job -- it's done */
  free(process_status[data_from]);
  process_status[data_from] = NULL;

  /* Now get the results */
  MPI_Recv (buf, bufsize, MPI_PACKED, data_from, SEARCH_STD_SCAN_RESULTS_TAG, 
	    MPI_COMM_WORLD, &status);
  MPI_Unpack (buf, bufsize, &position, &results_type, 1, 
	      MPI_CHAR, MPI_COMM_WORLD);

  if (results_type == SEARCH_STD_SCAN_RESULTS) {
    MPI_Unpack (buf, bufsize, &position, &num_results, 1, 
		MPI_INT, MPI_COMM_WORLD);
    if (num_results > 0 && cur_seq->results[(int)in_revcomp] == NULL) {
      cur_seq->results[(int)in_revcomp] = CreateResults (INIT_RESULTS);
    }
    for (i=0; i<num_results; i++) {
      MPI_Unpack(buf, bufsize, &position, &start, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack(buf, bufsize, &position, &stop, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack(buf, bufsize, &position, &bestr, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack(buf, bufsize, &position, &score,1, MPI_FLOAT, MPI_COMM_WORLD);
      /* Don't report hits from first D nucleotides in second overlapping seq
	 because it wasn't a full analysis for seqs ending there -- longer
	 seqs missed */
      if (index == 1 || stop > D)
	report_hit (index+start-1, index+stop-1, bestr, score, 
		    cur_seq->results[(int)in_revcomp]);
    }
    cur_seq->chunks_sent--;
  } else if (results_type == ALN_RESULTS) {
    ESL_ALLOC(tr, sizeof(Parsetree_t));
    /* Get size of the tree */
    MPI_Unpack (buf, bufsize, &position, &(tr->memblock), 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack (buf, bufsize, &position, &(tr->n), 1, MPI_INT, MPI_COMM_WORLD);
    /* Allocate it */
    ESL_ALLOC(tr->emitl, sizeof(int)*tr->n);
    ESL_ALLOC(tr->emitr, sizeof(int)*tr->n);
    ESL_ALLOC(tr->state, sizeof(int)*tr->n);
    ESL_ALLOC(tr->nxtl,  sizeof(int)*tr->n);
    ESL_ALLOC(tr->nxtr,  sizeof(int)*tr->n);
    ESL_ALLOC(tr->prv,   sizeof(int)*tr->n);
    ESL_ALLOC(tr->mode,  sizeof(int)*tr->n);
    ESL_ALLOC(tr->nalloc,sizeof(int)*tr->n);
    tr->nalloc = tr->n;

    /* Unpack it */
    MPI_Unpack (buf, bufsize, &position, tr->emitl, tr->n, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack (buf, bufsize, &position, tr->emitr, tr->n, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack (buf, bufsize, &position, tr->state, tr->n, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack (buf, bufsize, &position, tr->nxtl, tr->n, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack (buf, bufsize, &position, tr->nxtr, tr->n, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack (buf, bufsize, &position, tr->prv, tr->n, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack (buf, bufsize, &position, tr->mode, tr->n, MPI_INT, MPI_COMM_WORLD);
    cur_seq->results[(int)in_revcomp]->data[index].tr = tr;
    cur_seq->alignments_sent--;
  } else {
    Die ("Got result type %d when expecting SEARCH_STD_SCAN_RESULTS (%d) or ALN_RESULTS (%d)\n", results_type, SEARCH_STD_SCAN_RESULTS, ALN_RESULTS);
  }
  free(buf);
  return(cur_seq_index);
}
#endif


