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

#include "funcs.h"
#include "structs.h"

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

  if (MPI_Unpack(buf, n, pos, cm->ndtype,           nnodes,  MPI_CHAR, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  
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

  if (MPI_Pack(cm->ndtype,           nnodes,  MPI_CHAR, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");

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
  n += sz; 
  /* nodemap */

  if (MPI_Pack_size(nnodes,   MPI_CHAR, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += sz; 
  /* ndtype */

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
  n += sz; 
  /* nodemap */

  if (MPI_Pack_size(nnodes,   MPI_CHAR, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += sz; 
  /* ndtype */

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
  ESL_DPRINTF1(("cm_dsq_MPISend(): dsq has size %d\n", n));

  /* Make sure the buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) {
    void *tmp;
    ESL_RALLOC(*buf, tmp, sizeof(char) * n);
    *nalloc = n; 
  }
  ESL_DPRINTF1(("cm_dsq_MPISend(): buffer is ready\n"));

  /* Pack the status code, L and dsq into the buffer */
  position = 0;
  code     = (dsq == NULL) ? eslEOD : eslOK;
  if (MPI_Pack(&code, 1, MPI_INT, *buf, n, &position, comm) != 0) ESL_EXCEPTION(eslESYS, "mpi pack failed"); 
  if (dsq != NULL) {
    if (MPI_Pack(&L,  1,   MPI_INT,  *buf, n, &position, comm) != 0) ESL_EXCEPTION(eslESYS, "mpi pack failed"); 
    if (MPI_Pack(dsq, L+2, MPI_BYTE, *buf, n, &position, comm) != 0) ESL_EXCEPTION(eslESYS, "mpi pack failed"); 
  }
  ESL_DPRINTF1(("cm_dsq_MPISend(): dsq is packed into %d bytes\n", position));

  /* Send the packed profile to destination  */
  if (MPI_Send(*buf, n, MPI_PACKED, dest, tag, comm) != 0) ESL_EXCEPTION(eslESYS, "mpi send failed");
  ESL_DPRINTF1(("cm_dsq_MPISend(): dsq is sent.\n"));
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
  ESL_DPRINTF1(("cm_dsq_MPIRecv(): about to receive dsq.\n"));
  if (MPI_Recv(*buf, n, MPI_PACKED, source, tag, comm, &mpistatus) != 0) ESL_XEXCEPTION(eslESYS, "mpi recv failed");
  ESL_DPRINTF1(("cm_dsq_MPIRecv(): dsq has been received.\n"));

  /* Unpack it - where the first integer is a status code, OK or EOD */
  ESL_DPRINTF1(("cm_dsq_MPIRecv(): about to unpack dsq.\n"));
  pos = 0;
  if (MPI_Unpack       (*buf, n, &pos, &code,                   1, MPI_INT,           comm) != 0)  ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (code == eslEOD) { status = eslEOD; goto ERROR; }

  if (MPI_Unpack       (*buf, n, &pos, &L,                      1, MPI_INT,           comm) != 0)  ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  ESL_ALLOC(dsq, sizeof(ESL_DSQ) * (L+2));
  if (MPI_Unpack       (*buf, n, &pos, dsq,                  (L+2), MPI_BYTE,          comm) != 0)  ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  dsq[0] = dsq[(L+1)] = eslDSQ_SENTINEL; /* overwrite */
  ESL_DPRINTF1(("cm_dsq_MPIRecv(): dsq has been unpacked.\n"));

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
  ESL_DPRINTF1(("cm_search_results_MPISend(): results has size %d\n", n));

  /* Make sure the buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) {
    void *tmp;
    ESL_RALLOC(*buf, tmp, sizeof(char) * n);
    *nalloc = n; 
  }
  ESL_DPRINTF1(("cm_search_results_MPISend(): buffer is ready\n"));

  /* Pack the status code, and results into the buffer */
  position = 0;
  code     = (results == NULL) ? eslEOD : eslOK;
  if (MPI_Pack(&code, 1, MPI_INT, *buf, n, &position, comm) != 0) ESL_EXCEPTION(eslESYS, "mpi pack failed"); 
  if (results != NULL) {
    if ((status = cm_search_results_MPIPack(results, *buf, n, &position, comm)) != eslOK) return status;
  }
  ESL_DPRINTF1(("cm_search_results_MPISend(): results is packed into %d bytes\n", position));

  /* Send the packed profile to destination  */
  if (MPI_Send(*buf, n, MPI_PACKED, dest, tag, comm) != 0) ESL_EXCEPTION(eslESYS, "mpi send failed");
  ESL_DPRINTF1(("cm_search_results_MPISend(): results are sent.\n"));
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

  ESL_DPRINTF1(("cm_search_results_MPIPack(): ready.\n"));

  status = MPI_Pack((int *) &(results->num_results),   1, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  for (i = 0; i < results->num_results; i++) {
    status = cm_search_result_node_MPIPack(&(results->data[i]), buf, n, position, comm);  if (status != eslOK) return status;
  }
  ESL_DPRINTF1(("cm_search_results_MPIPack(): done. Packed %d bytes into buffer of size %d\n", *position, n));
  
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
  ESL_DPRINTF1(("cm_search_results_MPIUnpack(): %d results.\n", num_results));
  for (i = 0; i < num_results; i++) {
    status = cm_search_result_node_MPIUnpack(buf, n, pos, comm, &(results->data[i]));  if (status != eslOK) return status;
    assert(results->data[i].tr != NULL);
    assert(results->data[i].tr->n > 0);
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
  status = MPI_Pack_size (1, MPI_FLOAT, comm, &sz); n +=   sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  status = MPI_Pack_size (1, MPI_INT, comm, &sz);   n +=   sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  if(rnode->tr != NULL) {
    if ((status = cm_parsetree_MPIPackSize(rnode->tr, comm, &sz))  != eslOK) goto ERROR; n += sz;
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

  ESL_DPRINTF1(("cm_search_result_node_MPIPack(): ready.\n"));
  
  status = MPI_Pack((int *) &(rnode->start),   1, MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *) &(rnode->stop),    1, MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *) &(rnode->bestr),   1, MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *) &(rnode->score),   1, MPI_FLOAT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  if(rnode->tr != NULL) {
    has_tr = TRUE;
    status = MPI_Pack(&has_tr,                 1, MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
    status = cm_parsetree_MPIPack((Parsetree_t *) rnode->tr, buf, n, position, comm);  if (status != eslOK) return status;
  }
  else {
    has_tr = FALSE;
    status = MPI_Pack(&has_tr,                 1, MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  }

  ESL_DPRINTF1(("cm_search_result_node_MPIPack(): done. Packed %d bytes into buffer of size %d\n", *position, n));

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
  Parsetree_t *tr = NULL;

  status = MPI_Unpack (buf, n, pos, &start, 1, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &stop,  1, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &bestr, 1, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &score, 1, MPI_FLOAT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  rnode.start = start;
  rnode.stop  = stop;
  rnode.bestr = bestr;
  rnode.score = score;
  rnode.tr    = NULL;

  /* optionally, unpack a parsetree */
  status = MPI_Unpack (buf, n, pos, &has_tr, 1, MPI_INT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if(has_tr == TRUE) {
    status   = cm_parsetree_MPIUnpack(buf, n, pos, comm, &tr);  if (status != eslOK) return status;
    rnode.tr = tr;
  }
  
  *ret_rnode = rnode;
  return eslOK;
  
 ERROR:
  if(tr != NULL) FreeParsetree(tr);
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
 * Args:      seqs_to_aln - seqs_to_aln_t object to send
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
cm_seqs_to_aln_MPISend(seqs_to_aln_t *seqs_to_aln, int offset, int nseq_to_send, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc)
{
  int   status;
  int   code;
  int   sz, n, position;

  /* First, figure out the size of the work unit */
  if (MPI_Pack_size(2, MPI_INT, comm, &n) != 0) ESL_EXCEPTION(eslESYS, "mpi pack size failed"); 
  if (seqs_to_aln != NULL) { 
    if ((status = cm_seqs_to_aln_MPIPackSize(seqs_to_aln, offset, nseq_to_send, comm, &sz)) != eslOK) return status;
    n += sz;
  }
  ESL_DPRINTF1(("cm_seqs_to_aln_MPISend(): seqs_to_aln has size %d\n", n));

  /* Make sure the buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) {
    void *tmp;
    ESL_RALLOC(*buf, tmp, sizeof(char) * n);
    *nalloc = n; 
  }
  ESL_DPRINTF1(("cm_seqs_to_aln_MPISend(): buffer is ready\n"));

  /* Pack the status code, and seqs_to_aln into the buffer */
  position = 0;
  code     = (seqs_to_aln == NULL) ? eslEOD : eslOK;
  if (MPI_Pack(&code, 1, MPI_INT, *buf, n, &position, comm) != 0) ESL_EXCEPTION(eslESYS, "mpi pack failed"); 
  if (seqs_to_aln != NULL) {
    if ((status = cm_seqs_to_aln_MPIPack(seqs_to_aln, offset, nseq_to_send, *buf, n, &position, comm)) != eslOK) return status;
  }
  ESL_DPRINTF1(("cm_seqs_to_aln_MPISend(): seqs_to_aln is packed into %d bytes\n", position));

  /* Send the packed profile to destination  */
  if (MPI_Send(*buf, n, MPI_PACKED, dest, tag, comm) != 0) ESL_EXCEPTION(eslESYS, "mpi send failed");
  ESL_DPRINTF1(("cm_seqs_to_aln_MPISend(): seqs_to_aln are sent.\n"));
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
  ESL_DPRINTF1(("cm_seqs_to_aln_MPIRecv(): about to receive dsq.\n"));
  if (MPI_Recv(*buf, n, MPI_PACKED, source, tag, comm, &mpistatus) != 0) ESL_XEXCEPTION(eslESYS, "mpi recv failed");
  ESL_DPRINTF1(("cm_seqs_to_aln_MPIRecv(): dsq has been received.\n"));

  /* Unpack it - where the first integer is a status code, OK or EOD */
  ESL_DPRINTF1(("cm_seqs_to_aln_MPIRecv(): about to unpack dsq.\n"));
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

  status = MPI_Pack_size (1, MPI_INT, comm, &sz); n += 6*sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");

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

  ESL_DPRINTF1(("cm_seqs_to_aln_MPIPack(): ready.\n"));
  ESL_DASSERT1((seqs_to_aln->nseq >= (offset + nseq_to_pack)));

  /* determine what information we have in the seqs_to_aln object */
  int has_sq      = TRUE;
  int has_tr      = TRUE;
  int has_cp9_tr  = TRUE;
  int has_post    = TRUE;  
  int has_sc      = TRUE;  

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

  status = MPI_Pack((int *) &(nseq_to_pack), 1, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *) &(has_sq),       1, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *) &(has_tr),       1, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *) &(has_cp9_tr),   1, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *) &(has_post),     1, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *) &(has_sc),       1, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");

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
      /* we call PackOpt, even though we know we should have a valid postal code */
      status = esl_mpi_PackOpt(seqs_to_aln->postcode[i], -1, MPI_CHAR, buf, n, position,  comm); if (status != eslOK) return status;
    }

  if(has_sc)
    status = MPI_Pack(&(seqs_to_aln->sc[i]), nseq_to_pack, MPI_FLOAT, buf, n, position, comm); if (status != eslOK) return status;

  ESL_DPRINTF1(("cm_seqs_to_aln_MPIPack(): done. Packed %d bytes into buffer of size %d\n", *position, n));
  
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

  status = MPI_Unpack (buf, n, pos, &num_seqs_to_aln, 1, MPI_INT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &has_sq,          1, MPI_INT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &has_tr,          1, MPI_INT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &has_cp9_tr,      1, MPI_INT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &has_post,        1, MPI_INT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &has_sc,          1, MPI_INT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if(num_seqs_to_aln == 0) { status = eslOK; goto ERROR; }

  seqs_to_aln = CreateSeqsToAln(num_seqs_to_aln, FALSE);
  ESL_DPRINTF1(("cm_seqs_to_aln_MPIUnpack(): %d seqs_to_aln.\n", num_seqs_to_aln));

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

  ESL_DPRINTF1(("cm_parsetree_MPIPack(): ready.\n"));
  
  status = MPI_Pack((int *) &(tr->n),        1, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *) &(tr->memblock), 1, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(tr->emitl,           tr->n, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(tr->emitr,           tr->n, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(tr->state,           tr->n, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(tr->mode,            tr->n, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(tr->nxtl,            tr->n, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(tr->nxtr,            tr->n, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(tr->prv,             tr->n, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");

  ESL_DPRINTF1(("cm_parsetree_MPIPack(): done. Packed %d bytes into buffer of size %d\n", *position, n));

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

  ESL_DPRINTF1(("cm_cp9trace_MPIPack(): ready.\n"));
  
  status = MPI_Pack((int *) &(cp9_tr->tlen), 1, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(cp9_tr->statetype,       cp9_tr->tlen, MPI_CHAR, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(cp9_tr->nodeidx,         cp9_tr->tlen, MPI_INT,  buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(cp9_tr->pos,             cp9_tr->tlen, MPI_INT,  buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");

  ESL_DPRINTF1(("cm_cp9trace_MPIPack(): done. Packed %d bytes into buffer of size %d\n", *position, n));

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
  if(! (sq->flags  & eslSQ_DIGITAL)) ESL_XEXCEPTION(eslESYS, "cm_digitized_sq_MPIPackSize, sq not digitized.");

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
  if(! (sq->flags  & eslSQ_DIGITAL)) ESL_EXCEPTION(eslESYS, "cm_digitized_sq_MPIPackSize, sq not digitized.");

  ESL_DPRINTF1(("cm_digitized_sq_MPIPack(): ready.\n"));
  
  status = MPI_Pack((int *) &(sq->n),             1, MPI_INT, buf, n, position,  comm);  if (status != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  status = esl_mpi_PackOpt(sq->name, -1, MPI_CHAR, buf, n, position,  comm);             if (status != eslOK) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((ESL_DSQ *) sq->dsq, sq->n+2,MPI_BYTE, buf, n, position, comm);   if (status != 0)     ESL_EXCEPTION(eslESYS, "pack failed");

  ESL_DPRINTF1(("cm_digitized_sq_MPIPack(): done. Packed %d bytes into buffer of size %d\n", *position, n));

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
/*HEREHEREHERHEERERE old funcs below EPN, Tue Aug 28 15:51:13 2007*/
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
  ESL_DSQ *dsq = NULL;
  int   L;
  float best_sc;

  int doing_cm_stats  = FALSE;
  int doing_cp9_stats = FALSE;

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
				       NULL,  /* filter fraction, nobody cares */
				       FALSE);/* don't align hits */
      
      MPI_Send(&(best_sc), 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
      free(dsq);
    }
  return;

}

/* Function:  dsq_MPISend()
 * Incept:    EPN, Wed May  9 17:30:14 2007
 *
 * Purpose:   Send sequence <dsq> to processor <dest>.
 *            
 * Returns:   eslOK on success; eslEINVAL if <dsq> is NULL
 *            and eslESYS if there is an MPI error.
 */
int
dsq_MPISend(ESL_DSQ *dsq, int L, int dest)
{
  int status;

  if(dsq == NULL) { status = eslESYS; goto ERROR; }

  if(MPI_Send(&L, 1, MPI_INT, dest, 0, MPI_COMM_WORLD) != 0) ESL_EXCEPTION(eslESYS, "mpi send failed.");
  /* receiver will now allocate storage, before reading on...*/
  if(MPI_Send(dsq, (L+2), MPI_BYTE, dest, 0, MPI_COMM_WORLD) != 0) ESL_EXCEPTION(eslESYS, "mpi send failed.");
  return eslOK;

 ERROR: 
  return status;
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
dsq_MPIRecv(ESL_DSQ **ret_dsq, int *ret_L)
{
  int status;
  char *dsq = NULL;
  MPI_Status mpistatus;
  int L;
  
  if(MPI_Recv(&L, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpistatus) != 0) ESL_EXCEPTION(eslESYS, "mpi receive failed.");
  ESL_ALLOC(dsq, sizeof(ESL_DSQ) * (L+2));
  if(MPI_Recv(dsq, (L+2), MPI_CHAR, 0, 0, MPI_COMM_WORLD, &mpistatus) != 0); ESL_EXCEPTION(eslESYS, "mpi receive failed.");
  dsq[0] = dsq[(L+1)] = eslDSQ_SENTINEL;
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
  ESL_DSQ *dsq = NULL;
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
					   NULL,  /* filter fraction, nobody cares */
					   FALSE);/* don't align hits */
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
					     NULL,  /* filter fraction N/A */
					     FALSE);/* don't align hits */
	  
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
						 NULL,  /* filter fraction N/A */
						 FALSE);/* don't align hits */
	    }
	}
      else
	scores[0] = actually_search_target(cm, dsq, 1, L, 
					   0.,    /* minimum CM bit cutoff, irrelevant (?) */
					   0.,    /* minimum CP9 bit cutoff, irrelevant (?) */
					   NULL,  /* do not keep results */
					   FALSE, /* do not filter with a CP9 HMM */
					   FALSE, FALSE, /* not doing CM or CP9 Gumbel calcs */
					   NULL,  /* filter fraction, nobody cares */
					   FALSE);/* don't align hits */
      
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
	    SEARCH_STD_SEARCH_RESULTS_SIZE_TAG, MPI_COMM_WORLD, &status);
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
  MPI_Recv (buf, bufsize, MPI_PACKED, data_from, SEARCH_STD_SEARCH_RESULTS_TAG, 
	    MPI_COMM_WORLD, &status);
  MPI_Unpack (buf, bufsize, &position, &results_type, 1, 
	      MPI_CHAR, MPI_COMM_WORLD);

  if (results_type == SEARCH_STD_SEARCH_RESULTS) {
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
    Die ("Got result type %d when expecting SEARCH_STD_SEARCH_RESULTS (%d) or ALN_RESULTS (%d)\n", results_type, SEARCH_STD_SEARCH_RESULTS, ALN_RESULTS);
  }
  free(buf);
  return(cur_seq_index);
}
#endif


#if 0 
/* Unnecessary? maybe not, keep here temporarily */

/* Function:  cm_search_results_MPIRecv()
 *
 * Incept:    EPN, Wed Aug 29 05:05:06 2007
 *
 * Purpose:   Receives a work unit that consists of search results
 *            <source> (<0..nproc-1>, or <MPI_ANY_SOURCE>) tagged 
 *            as <tag> from communicator <comm>. 
 *            
 *            Work units are prefixed by a status code. If the unit's
 *            code is <eslOK> and no errors are encountered, this
 *            routine will return <eslOK> and a non-<NULL> <*ret_results>.
 *            If the unit's code is <eslEOD> (a shutdown signal), 
 *            this routine returns <eslEOD> and <*ret_results> is <NULL>.
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
 *            <eslEOD>, and <*ret_results> is <NULL>. <*ret_mpistatus> is filled
 *            with mpistatus.
 *            
 * Returns:   <eslOK> on success. <*ret_results> contains the new 
 *            results; they are allocated here, and the caller is 
 *            responsible for free'ing them.  <*buf> may have been 
 *            reallocated to a larger size, and <*nalloc> may have 
 *            been increased. <*ret_mpistatus> is filled with mpistatus.
 *
 *
 * Throws:    <eslESYS> if an MPI call fails; <eslEMEM> if an allocation fails.
 *            In either case, <*ret_results> is NULL, and the <buf> and its size
 *            <*nalloc> remain valid.
 */
int
cm_search_results_MPIRecv(int source, int tag, MPI_Comm comm, char **buf, int *nalloc, search_results_t  **ret_results, 
			  MPI_Status *ret_mpistatus)
{
  int         status, code;
  int         n;
  int         pos;

  /* Probe first, because we need to know if our buffer is big enough. */
  if (MPI_Probe(source, tag, comm, &mpistatus)  != 0) ESL_XEXCEPTION(eslESYS, "mpi probe failed");
  if (MPI_Get_count(ret_mpistatus, MPI_PACKED, &n) != 0) ESL_XEXCEPTION(eslESYS, "mpi get count failed");

  /* Make sure the buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) {
    void *tmp;
    ESL_RALLOC(*buf, tmp, sizeof(char) * n);
    *nalloc = n; 
  }

  /* Receive the packed work unit */
  ESL_DPRINTF1(("cm_search_results_MPIRecv(): about to receive results.\n"));
  if (MPI_Recv(*buf, n, MPI_PACKED, source, tag, comm, &mpistatus) != 0) ESL_XEXCEPTION(eslESYS, "mpi recv failed");
  ESL_DPRINTF1(("cm_search_results_MPIRecv(): results have been received.\n"));

  /* Unpack it - where the first integer is a status code, OK or EOD */
  ESL_DPRINTF1(("cm_search_results_MPIRecv(): about to unpack dsq.\n"));
  pos = 0;
  if (MPI_Unpack       (*buf, n, &pos, &code,                   1, MPI_INT,           comm) != 0)  ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (code == eslEOD) { status = eslEOD; goto ERROR; }

  return cm_search_results_MPIUnpack(*buf, *nalloc, &pos, comm, ret_results);

 ERROR:
  *ret_results = NULL;
  return status;
}
#endif
