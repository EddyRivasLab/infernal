/* mpisupport.c
 * Optional support for MPI parallelization in Infernal.
 * 
 * EPN, Mon Aug 27 12:38:13 2007
 * SVN $Id$
 */

#include "esl_config.h"
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
    /* if ((status = cm_MPIPackSize(cm, comm, &sz)) != eslOK) return status; */
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
    /*if ((status = cm_MPIPack(cm, *buf, n, &pos, comm)) != eslOK) return status;*/
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

  /*return cm_MPIUnpack(abc, *buf, *nalloc, &pos, comm, ret_cm);*/
  return cm_justread_MPIUnpack(abc, *buf, *nalloc, &pos, comm, ret_cm);

 ERROR:
  *ret_cm = NULL;
  return status;
}

/* Function:  cm_MPIUnpack()
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
  int has_stats;

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
  CMZero(cm);

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

  if (MPI_Unpack(buf, n, pos, &(cm->el_selfsc),          1, MPI_FLOAT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->enf_scdiff),         1, MPI_FLOAT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
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

  if (MPI_Unpack(buf, n, pos, &(has_stats),              1, MPI_INT,   comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");

  /* free the cm->comlog that was allocated inside the CreateCM() call, we'll allocate a new one, a bit messy */
  FreeComLog(cm->comlog);
  status   = comlog_MPIUnpack(buf, n, pos, comm, &(cm->comlog));  if (status != eslOK) return status;

  if (has_stats) { 
    status   = cmstats_MPIUnpack(buf, n, pos, comm, &(cm->stats));  if (status != eslOK) return status;
  }  

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

/* Function:  cm_justread_MPIUnpack()
 * Synopsis:  Unpacks a CM <cm> (packed by cm_justread_MPIPack()) from an MPI buffer.
 * Incept:    EPN, Mon Aug 27 14:44:11 2007
 *
 * Purpose:   Unpack a newly allocated CM from MPI packed buffer
 *            <buf>, starting from position <*pos>, where the total length
 *            of the buffer in bytes is <n>. Differs from cm_MPIPack()
 *            in that the full CM data structure was not packed, only
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
  int     status;
  CM_t *cm = NULL;
  int M, nnodes, K, atype;
  int has_stats;

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
  CMZero(cm);

  /* Unpack the rest of the CM */
  if (MPI_Unpack(buf, n, pos, &(cm->flags),              1,   MPI_INT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  /* note cm->flags is how it was immediately after CM is read from file, so the only flags that can possibly be raised are:
   * CMH_GA, CMH_TC, CMH_NC, CMH_GUMBEL_STATS and CMH_FILTER_STATS, which is okay b/c all that info is transmitted in this func */
  if (MPI_Unpack(buf, n, pos, &(cm->nseq),               1,   MPI_INT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->clen),               1,   MPI_INT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");

  if (MPI_Unpack(buf, n, pos, &(cm->el_selfsc),          1, MPI_FLOAT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->eff_nseq),           1, MPI_FLOAT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->ga),                 1, MPI_FLOAT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->tc),                 1, MPI_FLOAT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(cm->nc),                 1, MPI_FLOAT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, cm->null,                  K, MPI_FLOAT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");

  if (MPI_Unpack(buf, n, pos, cm->sttype,            (M+1),  MPI_CHAR, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, cm->stid,              (M+1),  MPI_CHAR, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");

  if (MPI_Unpack(buf, n, pos, cm->ndidx,                 M,   MPI_INT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, cm->cfirst,                M,   MPI_INT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, cm->cnum,                  M,   MPI_INT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, cm->plast,                 M,   MPI_INT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, cm->pnum,                  M,   MPI_INT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");

  if (MPI_Unpack(buf, n, pos, cm->nodemap,          nnodes,   MPI_INT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, cm->ndtype,           nnodes,  MPI_CHAR, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  
  if (MPI_Unpack(buf, n, pos, cm->e[0],              M*K*K, MPI_FLOAT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");

  if (MPI_Unpack(buf, n, pos, cm->t[0],       M*MAXCONNECT, MPI_FLOAT, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  
  /* free the cm->comlog that was allocated inside the CreateCM() call, we'll allocate a new one, a bit messy */
  FreeComLog(cm->comlog);
  status   = comlog_MPIUnpack(buf, n, pos, comm, &(cm->comlog));  if (status != eslOK) return status;

  if (MPI_Unpack(buf, n, pos, &(has_stats),              1, MPI_INT,   comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if (has_stats) { 
    status   = cmstats_MPIUnpack(buf, n, pos, comm, &(cm->stats));  if (status != eslOK) return status;
  }  

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
  cm_Fail("EPN, Fri Nov  9 08:55:23 2007, cm_MPIPack() shouldn't be used until oesc's, scanmatrix, and CM_HB_MX's are handled. Why not use cm_justread_MPIPack?()\n");

  int   status;
  int   K      = cm->abc->K;
  int   M      = cm->M;
  int   nnodes = cm->nodes;
  int   has_stats;

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

  if (MPI_Pack(&(cm->el_selfsc),          1, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&(cm->enf_scdiff),         1, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
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

  status = comlog_MPIPack(cm->comlog, buf, n, pos, comm);  if (status != eslOK) return status;

  has_stats = (cm->stats == NULL) ? FALSE : TRUE;
  if (MPI_Pack((int *) &(has_stats),      1, MPI_INT,   buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (has_stats) { 
    status = cmstats_MPIPack(cm->stats, buf, n, pos, comm);  if (status != eslOK) return status;
  }
  
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
  int   K      = cm->abc->K;
  int   M      = cm->M;
  int   nnodes = cm->nodes;
  int   has_stats;

  if (MPI_Pack(&M,                        1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(&nnodes,                   1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack((void *) &(cm->abc->type), 1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");

  if (MPI_Pack(&(cm->flags),              1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  /* note cm->flags is how it is immediately after CM is read from file, so the only flags that can possibly be raised are:
   * CMH_GA, CMH_TC, CMH_NC, CMH_GUMBEL_STATS and CMH_FILTER_STATS, which is okay b/c all the relevant info for those flags
   * is transmitted in this func */

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
  if (MPI_Pack(cm->ndtype,           nnodes,  MPI_CHAR, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");

  if (MPI_Pack(cm->e[0],              M*K*K, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");

  if (MPI_Pack(cm->t[0],       M*MAXCONNECT, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  
  status = comlog_MPIPack(cm->comlog, buf, n, pos, comm);  if (status != eslOK) return status;

  has_stats = (cm->stats == NULL) ? FALSE : TRUE;

  if (MPI_Pack((int *) &(has_stats),      1, MPI_INT,   buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (has_stats) { 
    status = cmstats_MPIPack(cm->stats, buf, n, pos, comm);  if (status != eslOK) return status;
  }

  if ((status = esl_mpi_PackOpt(cm->name,    -1, MPI_CHAR, buf, n, pos, comm)) != eslOK) return status; 
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
  n += 12*sz; 
  /* M, nodes, abc->type, flags, config_opts, search_opts, align_opts, nseq, clen, iel_selfsc, (10)
   * W, enf_start (2) 
   */

  if (MPI_Pack_size(1,       MPI_FLOAT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += (8+K)*sz; 
  /* el_selfsc, enf_scdiff, pbegin, pend, eff_nseq (5) 
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

  if ((status = comlog_MPIPackSize(cm->comlog, comm, &sz)) != eslOK) return status;
  n += sz;

  if (MPI_Pack_size(1, MPI_INT,comm,&sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /* has_stats: flag which tells worker whether or not it's about to get stats */

  if(cm->stats != NULL) { 
    if ((status = cmstats_MPIPackSize(cm->stats, comm, &sz)) != eslOK) return status;
    n += sz;
  }

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
  int   n = 0;
  int   K = cm->abc->K;
  int   M = cm->M;
  int   nnodes = cm->nodes;
  int   sz;

  if (MPI_Pack_size(1,         MPI_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");   
  n += 6*sz; /* M, nodes, abc->type, flags, nseq, clen */ 

  if (MPI_Pack_size(1,       MPI_FLOAT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += (5+K)*sz; 
  /* el_selfsc, eff_nseq, ga, tc, nc, null (5+K) */

  if (MPI_Pack_size((M+1),    MPI_CHAR, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += 2*sz; 
  /* sttype, stid */

  if (MPI_Pack_size(M,         MPI_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += 5*sz; 
  /* ndidx, cfirst, cnum, plast, pnum */

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

  if ((status = comlog_MPIPackSize(cm->comlog, comm, &sz)) != eslOK) return status;
  n += sz;

  if (MPI_Pack_size(1, MPI_INT,comm,&sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  n += sz;
  /* has_stats: flag which tells worker whether or not it's about to get stats */

  if(cm->stats != NULL) { 
    if ((status = cmstats_MPIPackSize(cm->stats, comm, &sz)) != eslOK) return status;
    n += sz;
  }

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
  if(rnode->pcode1 != NULL) { 
    if ((status = esl_mpi_PackOptSize(rnode->pcode1, -1, MPI_CHAR, comm, &sz)) != eslOK) goto ERROR; n += sz;
  }
  if(rnode->pcode2 != NULL) { 
    if ((status = esl_mpi_PackOptSize(rnode->pcode2, -1, MPI_CHAR, comm, &sz)) != eslOK) goto ERROR; n += sz;
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
  has_pcodes = (rnode->pcode1 != NULL && rnode->pcode2 != NULL) ? TRUE : FALSE;
  status = MPI_Pack(&has_tr,                 1, MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(&has_pcodes,             1, MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");

  if(has_tr) { 
    status = cm_parsetree_MPIPack((Parsetree_t *) rnode->tr, buf, n, position, comm);  if (status != eslOK) return status;
  }
  if(has_pcodes) {
    /* we call PackOpt, even though we know we should have valid postal codes */
    status = esl_mpi_PackOpt(rnode->pcode1, -1, MPI_CHAR, buf, n, position,  comm); if (status != eslOK) return status;
    status = esl_mpi_PackOpt(rnode->pcode2, -1, MPI_CHAR, buf, n, position,  comm); if (status != eslOK) return status;
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
  char *pcode1 = NULL;
  char *pcode2 = NULL;

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
  rnode.pcode1= NULL;
  rnode.pcode2= NULL;

  /* optionally, unpack a parsetree */
  if(has_tr) {
    status   = cm_parsetree_MPIUnpack(buf, n, pos, comm, &tr);  if (status != eslOK) return status;
    rnode.tr = tr;
  }
  if(has_pcodes) {
    status = esl_mpi_UnpackOpt(buf, n, pos, (void **) &(pcode1), NULL, MPI_CHAR, comm); if (status != eslOK) goto ERROR;;
    status = esl_mpi_UnpackOpt(buf, n, pos, (void **) &(pcode2), NULL, MPI_CHAR, comm); if (status != eslOK) goto ERROR;;
    rnode.pcode1 = pcode1;
    rnode.pcode2 = pcode2;
  }
  
  *ret_rnode = rnode;
  return eslOK;
  
 ERROR:
  if(tr != NULL) FreeParsetree(tr);
  if(pcode1 != NULL) free(pcode1);
  if(pcode2 != NULL) free(pcode2);
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

  if(seqs_to_aln->postcode1 == NULL) has_post = FALSE;
  else
    for (i = offset; i < offset + nseq_to_pack; i++) 
      if(seqs_to_aln->postcode1[i] == NULL) { has_post = FALSE; break; }

  if(seqs_to_aln->postcode2 == NULL) has_post = FALSE;
  else
    for (i = offset; i < offset + nseq_to_pack; i++) 
      if(seqs_to_aln->postcode2[i] == NULL) { has_post = FALSE; break; }

  if(seqs_to_aln->sc == NULL) has_sc = FALSE;
  else
    for (i = offset; i < offset + nseq_to_pack; i++) 
      if(! NOT_IMPOSSIBLE(seqs_to_aln->sc[i])) { has_sc = FALSE; break; }

  status = MPI_Pack_size (1, MPI_INT, comm, &sz); n += 6*sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /* nseq_to_pack, has_sq, has_tr, has_cp9_tr, has_post, has_sc */

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
      if ((status = esl_mpi_PackOptSize(seqs_to_aln->postcode1[i], -1, MPI_CHAR, comm, &sz)) != eslOK) goto ERROR; n += sz;
      if ((status = esl_mpi_PackOptSize(seqs_to_aln->postcode2[i], -1, MPI_CHAR, comm, &sz)) != eslOK) goto ERROR; n += sz;
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

  ESL_DPRINTF2(("cm_seqs_to_aln_MPIPack(): ready.\n"));
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

  if(seqs_to_aln->postcode1 == NULL) has_post = FALSE;
  else
    for (i = offset; i < offset + nseq_to_pack; i++) 
      if(seqs_to_aln->postcode1[i] == NULL) { has_post = FALSE; break; }

  if(seqs_to_aln->postcode2 == NULL) has_post = FALSE;
  else
    for (i = offset; i < offset + nseq_to_pack; i++) 
      if(seqs_to_aln->postcode2[i] == NULL) { has_post = FALSE; break; }

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
      /* we call PackOpt, even though we know we should have valid postal codes */
      status = esl_mpi_PackOpt(seqs_to_aln->postcode1[i], -1, MPI_CHAR, buf, n, position,  comm); if (status != eslOK) return status;
      status = esl_mpi_PackOpt(seqs_to_aln->postcode2[i], -1, MPI_CHAR, buf, n, position,  comm); if (status != eslOK) return status;
    }

  if(has_sc)
    status = MPI_Pack((seqs_to_aln->sc + (offset * sizeof(float))), nseq_to_pack, MPI_FLOAT, buf, n, position, comm); if (status != eslOK) return status;

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

  status = MPI_Unpack (buf, n, pos, &num_seqs_to_aln, 1, MPI_INT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &has_sq,          1, MPI_INT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &has_tr,          1, MPI_INT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &has_cp9_tr,      1, MPI_INT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &has_post,        1, MPI_INT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &has_sc,          1, MPI_INT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
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
    ESL_ALLOC(seqs_to_aln->postcode1, sizeof(char *) * num_seqs_to_aln);
    ESL_ALLOC(seqs_to_aln->postcode2, sizeof(char *) * num_seqs_to_aln);
    for (i = 0; i < num_seqs_to_aln; i++) {
      status = esl_mpi_UnpackOpt(buf, n, pos, (void **) &(seqs_to_aln->postcode1[i]), NULL, MPI_CHAR, comm); if (status != eslOK) goto ERROR;;
      status = esl_mpi_UnpackOpt(buf, n, pos, (void **) &(seqs_to_aln->postcode2[i]), NULL, MPI_CHAR, comm); if (status != eslOK) goto ERROR;;
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


/* Function:  cmstats_MPIPackSize()
 * Synopsis:  Calculates number of bytes needed to pack a 
 *            CMStats_t object. Follows 'Purpose' 
 *            of other *_MPIPackSize() functions above. 
 *
 * Incept:    EPN, Wed Dec 12 05:37:01 2007
 *
 * Returns:   <eslOK> on success, and <*ret_n> contains the answer.
 *
 * Throws:    <eslESYS> if an MPI call fails, and <*ret_n> is set to 0. 
 *
 * Note:      The sizing calls here need to stay matched up with
 *            the calls in <cmstats_MPIPack()>.
 */
int
cmstats_MPIPackSize(CMStats_t *cmstats, MPI_Comm comm, int *ret_n)
{
  int status;
  int sz;
  int n = 0;
  int np = cmstats->np;
  int g,p,f;
  assert(np > 0);

  status = MPI_Pack_size(1, MPI_INT, comm, &sz);  n += sz;     if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /* np */
  status = MPI_Pack_size(1, MPI_INT, comm, &sz);  n += np*2*sz;if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /* ps, pe arrays */
  status = MPI_Pack_size(1, MPI_INT, comm, &sz);  n += sz*GC_SEGMENTS; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /* gc2p int array */

  for(g = 0; g < GUM_NMODES; g++) { 
    for(p = 0; p < np; p++) { /* add size of Gumbel for each gum_mode/partition combo */
      if ((status = gumbel_info_MPIPackSize(cmstats->gumAA[g][p], comm, &sz))  != eslOK) goto ERROR; n += sz;
    }
  }   
  for(f = 0; f < FTHR_NMODES; f++) { /* add size of best filter info for each filter mode */
    if ((status = best_filter_info_MPIPackSize(cmstats->bfA[f], comm, &sz))  != eslOK) goto ERROR; n += sz;
  }

  *ret_n = n;
  return eslOK;

 ERROR:
  *ret_n = 0;
  return status;
}

/* Function:  cmstats_MPIPack()
 * Synopsis:  Packs CMStats_t <cmstats> into MPI buffer.
 *            See 'Purpose','Returns' and 'Throws'
 *            of other *_MPIPack()'s for more info.
 *
 * Incept:    EPN, Wed Dec 12 05:29:10 2007
 */
int
cmstats_MPIPack(CMStats_t *cmstats, char *buf, int n, int *position, MPI_Comm comm)
{
  int status;
  int np = cmstats->np;
  int gc, g, p, f;

  assert(np > 0);

  ESL_DPRINTF2(("cmstats_MPIPack(): ready.\n"));
  
  status = MPI_Pack((int *) &(cmstats->np),         1, MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(cmstats->ps,                   np, MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(cmstats->pe,                   np, MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  for(gc = 0; gc < GC_SEGMENTS; gc++) { /* there must be a better way to do this, this is safe, slow route */
    status = MPI_Pack((int *) &(cmstats->gc2p[gc]), 1, MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  }
  for(g = 0; g < GUM_NMODES; g++) { 
    for(p = 0; p < np; p++) { /* pack gumbel for this gum_mode/partition combo */
      status = gumbel_info_MPIPack(cmstats->gumAA[g][p], buf, n, position, comm);  if (status != eslOK) return status;
    }
  }
  for(f = 0; f < FTHR_NMODES; f++) { 
    status = best_filter_info_MPIPack(cmstats->bfA[f], buf, n, position, comm);  if (status != eslOK) return status;
  }

  if (*position > n) ESL_EXCEPTION(eslEMEM, "buffer overflow");
  return eslOK;
}

/* Function:  cmstats_MPIUnpack()
 * Synopsis:  Unpacks CMStats_t <cmstats> from an MPI buffer.
 *            Follows 'Purpose', 'Returns', 'Throws' of other
 *            *_MPIUnpack() functions above.
 *
 * Incept:    EPN, Wed Dec 12 05:29:05 2007
 */
int
cmstats_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, CMStats_t **ret_cmstats)
{
  int status;
  CMStats_t *cmstats;
  int gc, g, f, p;
  int np;

  ESL_ALLOC(cmstats, sizeof(CMStats_t));
  status = MPI_Unpack (buf, n, pos, &(cmstats->np), 1, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  np = cmstats ->np;
  assert(np > 0);

  ESL_ALLOC(cmstats->ps, sizeof(int) * np);
  ESL_ALLOC(cmstats->pe, sizeof(int) * np);

  status = MPI_Unpack (buf, n, pos, cmstats->ps, np, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, cmstats->pe, np, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  for(gc = 0; gc < GC_SEGMENTS; gc++) { /* there must be a better way to do this, this is safe, slow route */
    status = MPI_Unpack (buf, n, pos, &(cmstats->gc2p[gc]), 1, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  }

  ESL_ALLOC(cmstats->gumAA, sizeof(GumbelInfo_t **) * GUM_NMODES);
  for(g = 0; g < GUM_NMODES; g++) { 
    ESL_ALLOC(cmstats->gumAA[g], sizeof(GumbelInfo_t *) * np);
    for(p = 0; p < np; p++) { /* pack gumbel for this gum_mode/partition combo */
      status = gumbel_info_MPIUnpack(buf, n, pos, comm, &(cmstats->gumAA[g][p]));  if (status != eslOK) return status;
    }
  }
  ESL_ALLOC(cmstats->bfA, sizeof(BestFilterInfo_t *) * FTHR_NMODES);
  for(f = 0; f < FTHR_NMODES; f++) { 
    status = best_filter_info_MPIUnpack(buf, n, pos, comm, &(cmstats->bfA[f]));  if (status != eslOK) return status;
  }
  *ret_cmstats = cmstats;
  return eslOK;

 ERROR:
  if(cmstats != NULL) FreeCMStats(cmstats);
  *ret_cmstats = NULL;
  return status;
}

/* Function:  gumbel_info_MPIPackSize()
 * Synopsis:  Calculates number of bytes needed to pack a 
 *            GumbelInfo_t object. Follows 'Purpose' 
 *            of other *_MPIPackSize() functions above. 
 *
 * Incept:    EPN, Wed Dec 12 05:00:01 2007
 *
 * Returns:   <eslOK> on success, and <*ret_n> contains the answer.
 *
 * Throws:    <eslESYS> if an MPI call fails, and <*ret_n> is set to 0. 
 *
 * Note:      The sizing calls here need to stay matched up with
 *            the calls in <gumbel_info_MPIPack()>.
 */
int
gumbel_info_MPIPackSize(GumbelInfo_t *gum, MPI_Comm comm, int *ret_n)
{
  int status;
  int sz;
  int n = 0;

  status = MPI_Pack_size(1, MPI_INT, comm, &sz);    n += 3*sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /* N, L, is_valid */
  status = MPI_Pack_size(1, MPI_DOUBLE, comm, &sz); n += 2*sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /* mu, lambda */

  *ret_n = n;
  return eslOK;

 ERROR:
  *ret_n = 0;
  return status;
}

/* Function:  gumbel_info_MPIPack()
 * Synopsis:  Packs GumbelInfo_t <gum> into MPI buffer.
 *            See 'Purpose','Returns' and 'Throws'
 *            of other *_MPIPack()'s for more info.
 *
 * Incept:    EPN, Wed Dec 12 05:03:40 2007
 */
int
gumbel_info_MPIPack(GumbelInfo_t *gum, char *buf, int n, int *position, MPI_Comm comm)
{
  int status;

  ESL_DPRINTF2(("gumbel_info_MPIPack(): ready.\n"));
  
  status = MPI_Pack((int *) &(gum->N),        1, MPI_INT,    buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *) &(gum->L),        1, MPI_INT,    buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((double *) &(gum->mu),    1, MPI_DOUBLE, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((double *) &(gum->lambda),1, MPI_DOUBLE, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *) &(gum->is_valid), 1, MPI_INT,    buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  ESL_DPRINTF2(("gumbel_info_results_MPIPack(): done. Packed %d bytes into buffer of size %d\n", *position, n));

  if (*position > n) ESL_EXCEPTION(eslEMEM, "buffer overflow");
  return eslOK;
}

/* Function:  gumbel_info_MPIUnpack()
 * Synopsis:  Unpacks GumbelInfo_t<gum> from an MPI buffer.
 *            Follows 'Purpose', 'Returns', 'Throws' of other
 *            *_MPIUnpack() functions above.
 *
 * Incept:    EPN, Wed Dec 12 05:29:19 2007
 */
int
gumbel_info_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, GumbelInfo_t **ret_gum)
{
  int status;
  GumbelInfo_t *gum;

  ESL_ALLOC(gum, sizeof(GumbelInfo_t));
  status = MPI_Unpack (buf, n, pos, &(gum->N),        1, MPI_INT,    comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &(gum->L),        1, MPI_INT,    comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &(gum->mu),       1, MPI_DOUBLE, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &(gum->lambda),   1, MPI_DOUBLE, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &(gum->is_valid), 1, MPI_INT,    comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  *ret_gum = gum;
  return eslOK;

 ERROR:
  if(gum != NULL) free(gum);
  *ret_gum = NULL;
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
  int is_hybrid;
  int p;

  is_hybrid = (bf->ftype == FILTER_WITH_HYBRID) ? TRUE : FALSE;
  if(!is_hybrid) assert(bf->np == 0);

  status = MPI_Pack_size(1, MPI_INT, comm, &sz);    n += 6*sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /* cm_M, N, db_size, is_valid, ftype, np */
  status = MPI_Pack_size(1, MPI_FLOAT, comm, &sz);  n += 6*sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /* cm_eval, F, e_cutoff, full_cm_ncalcs, fil_ncalcs, fil_plus_surv_ncalcs */
  if(is_hybrid) { /* v_isroot and hgumA are valid, we pack them */
    status = MPI_Pack_size(1, MPI_INT, comm, &sz);  n += bf->np*sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
    /* v_isroot array */
    for(p = 0; p < bf->np; p++) { /* add size of Gumbel for each partition */
      if ((status = gumbel_info_MPIPackSize(bf->hgumA[p], comm, &sz))  != eslOK) goto ERROR; n += sz;
    }
  }    
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
  int p;

  ESL_DPRINTF2(("best_filter_info_MPIPack(): ready.\n"));
  
  status = MPI_Pack((int *) &(bf->cm_M),                  1, MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *) &(bf->N),                     1, MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *) &(bf->db_size),               1, MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *) &(bf->is_valid),              1, MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *) &(bf->ftype),                 1, MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((int *) &(bf->np),                    1, MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");

  status = MPI_Pack((float *) &(bf->cm_eval),             1, MPI_FLOAT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((float *) &(bf->F),                   1, MPI_FLOAT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((float *) &(bf->e_cutoff),            1, MPI_FLOAT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((float *) &(bf->full_cm_ncalcs),      1, MPI_FLOAT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((float *) &(bf->fil_ncalcs),          1, MPI_FLOAT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack((float *) &(bf->fil_plus_surv_ncalcs),1, MPI_FLOAT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");

  int is_hybrid = (bf->ftype == FILTER_WITH_HYBRID) ? TRUE : FALSE;
  if(is_hybrid) { /* get v_isroot and hgumA */
    assert((bf->np > 0));
    status = MPI_Pack(bf->v_isroot, bf->np, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
    for (p = 0; p < bf->np; p++) {
      status = gumbel_info_MPIPack(bf->hgumA[p], buf, n, position, comm);  if (status != eslOK) return status;
    }
  }

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
  int p;

  ESL_ALLOC(bf, sizeof(BestFilterInfo_t));
  status = MPI_Unpack (buf, n, pos, &(bf->cm_M),                1, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &(bf->N),                   1, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &(bf->db_size),             1, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &(bf->is_valid),            1, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &(bf->ftype),               1, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &(bf->np),                  1, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  status = MPI_Unpack (buf, n, pos, &(bf->cm_eval),             1, MPI_FLOAT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &(bf->F),                   1, MPI_FLOAT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &(bf->e_cutoff),            1, MPI_FLOAT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &(bf->full_cm_ncalcs),      1, MPI_FLOAT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &(bf->fil_ncalcs),          1, MPI_FLOAT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = MPI_Unpack (buf, n, pos, &(bf->fil_plus_surv_ncalcs),1, MPI_FLOAT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  int is_hybrid = (bf->ftype == FILTER_WITH_HYBRID) ? TRUE : FALSE;
  if(is_hybrid) { /* unpack v_isroot, hgumA */
    assert(bf->np > 0);
    ESL_ALLOC(bf->v_isroot, sizeof(int) * bf->np);
    ESL_ALLOC(bf->hgumA,    sizeof(GumbelInfo_t *) * bf->np);
    status = MPI_Unpack (buf, n, pos, bf->v_isroot,     bf->np, MPI_INT, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
    for(p = 0; p < bf->np; p++) { 
      status = gumbel_info_MPIUnpack(buf, n, pos, comm, &(bf->hgumA[p]));  if (status != eslOK) return status;
    }
  }    
  else { /* not a hybrid, point v_isroot, hgumA to NULL */
    bf->v_isroot = NULL;
    bf->hgumA = NULL;
  }

  *ret_bf = bf;
  return eslOK;

 ERROR:
  if(bf != NULL) free(bf);
  *ret_bf = NULL;
  return status;
}


/* Function:  cmcalibrate_cm_gumbel_results_MPIPackSize()
 * Synopsis:  Calculates number of bytes needed to pack a 
 *            results for CM scan for cmcalibrate.
 * Incept:    EPN, Thu Dec  6 16:44:17 2007
 *
 * Purpose:   Calculate an upper bound on the number of bytes
 *            that <cmcalibrate_cm_gumbel_results_MPIPack()> will need 
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
 *            the calls in <cmcalibrate_cm_gumbel_results_MPIPack()>.
 */
int
cmcalibrate_cm_gumbel_results_MPIPackSize(float **vscAA, int nseq, int M, MPI_Comm comm, int *ret_n)
{
  int status;
  int sz;
  int n = 0;

  status = MPI_Pack_size(1, MPI_INT, comm, &sz);   n += sz;      if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  status = MPI_Pack_size(M, MPI_FLOAT, comm, &sz); n += nseq*sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");

  *ret_n = n;
  return eslOK;

 ERROR:
  *ret_n = 0;
  return status;
}

/* Function:  cmcalibrate_cm_gumbel_results_MPIPack()
 * Synopsis:  Packs CM vscAA scores into MPI buffer.
 * Incept:    EPN, Thu Dec  6 16:47:58 2007
 *
 * Purpose:   Packs <vscAA> into an MPI packed message 
 *            buffer <buf> of length <n> bytes, 
 *            starting at byte position
 *            <*position>, for MPI communicator <comm>.
 *
 *            Note: <vscAA> is a 2D array, vscAA[0..v..M-1][0..i..nseq-1]
 *            holding the best score for each subtree rooted 
 *            at v for a CM scan (CYK/Inside) of sequence i.
 *            But we send it as a 1D array, vscA, of M * nseq floats,
 *            0..nseq-1 correspond to v==0, nseq..(2*nseq-1) correspond
 *            to v==1, etc.
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
cmcalibrate_cm_gumbel_results_MPIPack(float **vscAA, int nseq, int M, char *buf, int n, int *position, MPI_Comm comm)
{
  int status;
  int i,v,idx;
  float *vscA = NULL;

  ESL_DPRINTF2(("cmcalibrate_cm_gumbel_results_MPIPack(): ready.\n"));

  ESL_ALLOC(vscA, sizeof(float) * (M*nseq));
  idx = 0;
  for(v = 0; v < M; v++) 
    for(i = 0; i < nseq; i++)
      vscA[idx++] = vscAA[v][i];

  status = MPI_Pack((int *) &(nseq), 1, MPI_INT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(vscA,  (M*nseq), MPI_FLOAT,  buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");

  ESL_DPRINTF2(("cmcalibrate_cm_gumbel_results_MPIPack(): done. Packed %d bytes into buffer of size %d\n", *position, n));

  if (*position > n) ESL_EXCEPTION(eslEMEM, "buffer overflow");
  return eslOK;

 ERROR:
  if(vscA  != NULL) free(vscA);
  return status;
}

/* Function:  cmcalibrate_cm_gumbel_results_MPIUnpack()
 * Synopsis:  Unpacks <vscAA> from an MPI buffer.
 * Incept:    EPN, Wed Aug 29 05:10:20 2007
 *
 * Purpose:   Unpack a newly allocated set of scores <vscAA> from MPI packed buffer
 *            <buf>, starting from position <*pos>, where the total length
 *            of the buffer in bytes is <n>. 
 *
 *            Note: We return <ret_vscAA> as a 2D array, 
 *            ret_vscAA[0..v..M-1][0..i..nseq-1]
 *            holding the best score for each subtree rooted 
 *            at v for a CM scan (CYK/Inside) of sequence i.
 *            But vscA is sent as a 1D array, of M * nseq floats,
 *            0..nseq-1 correspond to v==0, nseq..(2*nseq-1) correspond
 *            to v==1, etc.
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
cmcalibrate_cm_gumbel_results_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, int M, float ***ret_vscAA, int *ret_nseq)
{
  int status;
  float  *vscA  = NULL;
  float **vscAA = NULL;
  int nseq = 0;
  int i, v, idx;

  status = MPI_Unpack (buf, n, pos, &nseq,        1, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  ESL_ALLOC(vscA, sizeof(float) * (M*nseq));
  status = MPI_Unpack (buf, n, pos, vscA, (M*nseq), MPI_FLOAT,  comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  ESL_ALLOC(vscAA, sizeof(float *) * (M));
  idx = 0;
  for(v = 0; v < M; v++) { 
    ESL_ALLOC(vscAA[v], sizeof(float) * nseq);
    for(i = 0; i < nseq; i++)
      vscAA[v][i] = vscA[idx++];
  }
  ESL_DASSERT1((idx == (M*nseq)));

  free(vscA);
  *ret_vscAA = vscAA;
  *ret_nseq = nseq;
  return eslOK;

 ERROR:
  if(vscA  != NULL) free(vscA);
  if(vscAA != NULL) { 
    for(i = 0; i < nseq; i++)
      free(vscAA[i]);
    free(vscAA);
  }
  *ret_vscAA = NULL;
  *ret_nseq = 0;
  return status;
}


/* Function:  cmcalibrate_cp9_gumbel_results_MPIPackSize()
 * Synopsis:  Calculates number of bytes needed to pack a 
 *            results for CM scan for cmcalibrate.
 * Incept:    EPN, Thu Dec  6 16:56:27 2007
 *
 * Purpose:   Calculate an upper bound on the number of bytes
 *            that <cmcalibrate_cp9_gumbel_results_MPIPack()> will need 
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
 *            the calls in <cmcalibrate_cp9_gumbel_results_MPIPack()>.
 */
int
cmcalibrate_cp9_gumbel_results_MPIPackSize(float *cp9scA, int nseq, MPI_Comm comm, int *ret_n)
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

/* Function:  cmcalibrate_cp9_gumbel_results_MPIPack()
 * Synopsis:  Packs CM vscAA scores into MPI buffer.
 * Incept:    EPN, Thu Dec  6 16:56:31 2007
 *
 * Purpose:   Packs <vscAA> into an MPI packed message 
 *            buffer <buf> of length <n> bytes, 
 *            starting at byte position
 *            <*position>, for MPI communicator <comm>.
 *
 *            <cp9scA> is an array, cp9scA[0..i..nseq-1]
 *            holding the best score for a CP9 scan (Viterbi
 *            or Forward) against sequence i.
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
cmcalibrate_cp9_gumbel_results_MPIPack(float *cp9scA, int nseq, char *buf, int n, int *position, MPI_Comm comm)
{
  int status;

  ESL_DPRINTF2(("cmcalibrate_cp9_gumbel_results_MPIPack(): ready.\n"));
  
  status = MPI_Pack((int *) &(nseq), 1, MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(cp9scA,       nseq, MPI_FLOAT, buf, n, position,  comm);     if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  ESL_DPRINTF2(("cmcalibrate_cp9_gumbel_results_MPIPack(): done. Packed %d bytes into buffer of size %d\n", *position, n));

  if (*position > n) ESL_EXCEPTION(eslEMEM, "buffer overflow");
  return eslOK;
}

/* Function:  cmcalibrate_cp9_gumbel_results_MPIUnpack()
 * Synopsis:  Unpacks <vscAA> from an MPI buffer.
 * Incept:    EPN, Thu Dec  6 16:56:36 2007
 *
 * Purpose:   Unpack a newly allocated set of cp9 scores <cp9scA> from MPI packed buffer
 *            <buf>, starting from position <*pos>, where the total length
 *            of the buffer in bytes is <n>. 
 *
 *            <cp9scA> is an array, cp9scA[0..i..nseq-1]
 *            holding the best score for a CP9 scan (Viterbi
 *            or Forward) against sequence i.
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
cmcalibrate_cp9_gumbel_results_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, float **ret_cp9scA, int *ret_nseq)
{
  int status;
  float *cp9scA;
  int nseq = 0;

  status = MPI_Unpack (buf, n, pos, &nseq,        1, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  ESL_ALLOC(cp9scA, sizeof(float) * nseq);
  status = MPI_Unpack (buf, n, pos, cp9scA, nseq, MPI_FLOAT,  comm);  if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  *ret_cp9scA = cp9scA;
  *ret_nseq   = nseq;
  return eslOK;

 ERROR:
  if(cp9scA != NULL) free(cp9scA);
  *ret_cp9scA = NULL;
  *ret_nseq = 0;
  return status;
}

/* Function:  cmcalibrate_cp9_filter_results_MPIPackSize()
 * Synopsis:  Calculates number of bytes needed to pack 
 *            CP9 filter results for cmcalibrate without --hybrid.
 *            enabled. 
 *           
 *            Differs from cmcalibrate_cp9_filter_results_hyb_MPIPackSize()
 *            in that vscAA, best scores for each state from CM scans 
 *            ARE NOT packed. They're irrelevant unless --hybrid was enabled.
 *
 *            Follows, 'Purpose', 'Returns', 'Throws' of
 *            the many other *_MPIPackSize() funcs above.
 *            
 * Incept:    EPN, Tue Jan  8 15:14:17 2008
 *           
 * Note:      The sizing calls here need to stay matched up with
 *            the calls in <cmcalibrate_cp9_filter_results_MPIPack()>.
 */
int
cmcalibrate_cp9_filter_results_MPIPackSize(int nseq, MPI_Comm comm, int *ret_n)
{
  int status;
  int sz;
  int n = 0;

  status = MPI_Pack_size(1, MPI_INT, comm, &sz);        n += sz;   if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /* for nseq */
  status = MPI_Pack_size(nseq, MPI_FLOAT, comm, &sz);   n += 2*sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /* for vit_cp9scA, fwd_cp9scA */
  status = MPI_Pack_size(nseq, MPI_INT,   comm, &sz);   n += sz;   if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /* for partA */

  *ret_n = n;
  return eslOK;

 ERROR:
  *ret_n = 0;
  return status;
}

/* Function:  cmcalibrate_cp9_filter_results_MPIPack()
 * Synopsis:  Packs cmcalibrate CP9 filter results into MPI buffer.
 *            Follows, 'Purpose', 'Returns', 'Throws' of
 *            the many other *_MPIPack() funcs above.
 * 
 * Incept:    EPN, Wed Dec 12 16:36:02 2007
 *
 *            Differs from cmcalibrate_cp9_filter_results_hyb_MPIPackSize()
 *            in that vscAA, best scores for each state from CM scans 
 *            ARE NOT packed, those scores are only relevant if --hybrid
 *            was enabled in cmcalibrate.
 *
 */
int
cmcalibrate_cp9_filter_results_MPIPack(float *vit_cp9scA, float *fwd_cp9scA, int *partA, int nseq, char *buf, int n, int *position, MPI_Comm comm)
{
  int status;

  ESL_DPRINTF2(("cmcalibrate_cp9_filter_results_MPIPack(): ready.\n"));

  status = MPI_Pack((int *) &(nseq), 1,        MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(vit_cp9scA,      nseq,     MPI_FLOAT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(fwd_cp9scA,      nseq,     MPI_FLOAT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(partA,           nseq,     MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");

  ESL_DPRINTF2(("cmcalibrate_cp9_filter_results_MPIPack(): done. Packed %d bytes into buffer of size %d\n", *position, n));

  if (*position > n) ESL_EXCEPTION(eslEMEM, "buffer overflow");
  return eslOK;

}

/* Function:  cmcalibrate_cp9_filter_results_MPIUnpack()
 * Synopsis:  Unpacks cmcalibrate cp9 filter results from an MPI buffer.
 *            Follows, 'Purpose', 'Returns', 'Throws' of
 *            the many other *_MPIUnpack() funcs above.
 * Incept:    EPN, Wed Dec 12 16:38:15 2007
 *
 *            Differs from cmcalibrate_cp9_filter_results_hyb_MPIPackSize()
 *            in that vscAA, best scores for each state from CM scans 
 *            is NOT packed. Those scores are only relevant if --hybrid
 *            was enabled for cmcalibrate.
 *
 */
int
cmcalibrate_cp9_filter_results_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, float **ret_vit_cp9scA, float **ret_fwd_cp9scA, int **ret_partA, int *ret_nseq)
{
  int status;
  float  *vit_cp9scA  = NULL;
  float  *fwd_cp9scA  = NULL;
  int    *partA       = NULL;
  int nseq = 0;

  status = MPI_Unpack (buf, n, pos, &nseq,        1, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  ESL_ALLOC(vit_cp9scA, sizeof(float) * nseq);
  status = MPI_Unpack (buf, n, pos, vit_cp9scA, nseq, MPI_FLOAT,  comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  ESL_ALLOC(fwd_cp9scA, sizeof(float) * nseq);
  status = MPI_Unpack (buf, n, pos, fwd_cp9scA, nseq, MPI_FLOAT,  comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  ESL_ALLOC(partA,      sizeof(int) * nseq);
  status = MPI_Unpack (buf, n, pos, partA, nseq, MPI_INT,    comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  *ret_vit_cp9scA = vit_cp9scA;
  *ret_fwd_cp9scA = fwd_cp9scA;
  *ret_partA      = partA;
  *ret_nseq = nseq;
  return eslOK;

  ESL_DPRINTF1(("cmcalibrate_cp9_filter_results_MPIUnpack() done.\n"));

 ERROR:
  if(vit_cp9scA != NULL) free(vit_cp9scA);
  if(fwd_cp9scA != NULL) free(fwd_cp9scA);
  if(partA      != NULL) free(partA);
  *ret_nseq = 0;
  return status;
}

/* Function:  cmcalibrate_cp9_filter_results_hyb_MPIPackSize()
 * Synopsis:  Calculates number of bytes needed to pack 
 *            CP9 filter results for cmcalibrate with --hybrid.
 *            enabled. 
 *           
 *            Differs from cmcalibrate_cp9_filter_results_MPIPackSize()
 *            in that vscAA, best scores for each state from CM scans 
 *            is packed. these scores are eventually used to calculate
 *            Gumbels for each possible sub CM state.
 *
 *            Follows, 'Purpose', 'Returns', 'Throws' of
 *            the many other *_MPIPackSize() funcs above.
 *            
 * Incept:    EPN, Wed Dec 12 16:30:20 2007
 *           
 * Note:      The sizing calls here need to stay matched up with
 *            the calls in <cmcalibrate_cp9_filter_results_hyb_MPIPack()>.
 */
int
cmcalibrate_cp9_filter_results_hyb_MPIPackSize(int nseq, int M, MPI_Comm comm, int *ret_n)
{
  int status;
  int sz;
  int n = 0;

  status = MPI_Pack_size(1, MPI_INT, comm, &sz);        n += sz;   if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /* for nseq */
  status = MPI_Pack_size(M*nseq, MPI_FLOAT, comm, &sz); n += sz;   if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /* for vscAA (we'll send it as a 1D array of M * nseq floats */
  status = MPI_Pack_size(nseq, MPI_FLOAT, comm, &sz);   n += 2*sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /* for vit_cp9scA, fwd_cp9scA */
  status = MPI_Pack_size(nseq, MPI_INT,   comm, &sz);   n += sz;   if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  /* for partA */

  *ret_n = n;
  return eslOK;

 ERROR:
  *ret_n = 0;
  return status;
}

/* Function:  cmcalibrate_cp9_filter_results_hyb_MPIPack()
 * Synopsis:  Packs cmcalibrate CP9 filter results into MPI buffer.
 *            Follows, 'Purpose', 'Returns', 'Throws' of
 *            the many other *_MPIPack() funcs above.
 * 
 * Incept:    EPN, Wed Dec 12 16:36:02 2007
 *
 *            Differs from cmcalibrate_cp9_filter_results_MPIPackSize()
 *            in that vscAA, best scores for each state from CM scans 
 *            is packed. these scores are eventually used to calculate
 *            Gumbels for each possible sub CM state.
 *
 *            Note: <vscAA> is a 2D array, vscAA[0..v..M-1][0..i..nseq-1]
 *            holding the best score for each subtree rooted 
 *            at v for a CM scan of sequence i.
 *            But we send it as a 1D array, vscA, of M * nseq floats,
 *            0..nseq-1 correspond to v==0, nseq..(2*nseq-1) correspond
 *            to v==1, etc.
 *
 */
int
cmcalibrate_cp9_filter_results_hyb_MPIPack(float **vscAA, float *vit_cp9scA, float *fwd_cp9scA, int *partA, int nseq, int M, char *buf, int n, int *position, MPI_Comm comm)
{
  int status;
  int i,v,idx;
  float *vscA = NULL;

  ESL_DPRINTF2(("cmcalibrate_cp9_filter_results_hyb_MPIPack(): ready.\n"));

  ESL_ALLOC(vscA, sizeof(float) * (M*nseq));
  idx = 0;
  for(v = 0; v < M; v++) 
    for(i = 0; i < nseq; i++)
      vscA[idx++] = vscAA[v][i];

  status = MPI_Pack((int *) &(nseq), 1,        MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(vscA,            (M*nseq), MPI_FLOAT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(vit_cp9scA,      nseq,     MPI_FLOAT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(fwd_cp9scA,      nseq,     MPI_FLOAT, buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");
  status = MPI_Pack(partA,           nseq,     MPI_INT,   buf, n, position,  comm); if (status != 0) ESL_EXCEPTION(eslESYS, "pack failed");

  ESL_DPRINTF2(("cmcalibrate_cp9_filter_results_hyb_MPIPack(): done. Packed %d bytes into buffer of size %d\n", *position, n));

  if (*position > n) ESL_EXCEPTION(eslEMEM, "buffer overflow");
  return eslOK;

 ERROR:
  if(vscA  != NULL) free(vscA);
  return status;
}

/* Function:  cmcalibrate_cp9_filter_results_hyb_MPIUnpack()
 * Synopsis:  Unpacks cmcalibrate cp9 filter results from an MPI buffer.
 *            Follows, 'Purpose', 'Returns', 'Throws' of
 *            the many other *_MPIUnpack() funcs above.
 * Incept:    EPN, Wed Dec 12 16:38:15 2007
 *
 *            Differs from cmcalibrate_cp9_filter_results_MPIPackSize()
 *            in that vscAA, best scores for each state from CM scans 
 *            is packed. these scores are eventually used to calculate
 *            Gumbels for each possible sub CM state.
 *
 *            Note: We return <ret_vscAA> as a 2D array, 
 *            ret_vscAA[0..v..M-1][0..i..nseq-1]
 *            holding the best score for each subtree rooted 
 *            at v for a CM scan of sequence i.
 *            But vscA is sent as a 1D array, of M * nseq floats,
 *            0..nseq-1 correspond to v==0, nseq..(2*nseq-1) correspond
 *            to v==1, etc.
 *
 */
int
cmcalibrate_cp9_filter_results_hyb_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, int M, float ***ret_vscAA, float **ret_vit_cp9scA, float **ret_fwd_cp9scA, int **ret_partA, int *ret_nseq)
{
  int status;
  float  *vit_cp9scA  = NULL;
  float  *fwd_cp9scA  = NULL;
  int    *partA       = NULL;
  float  *vscA  = NULL;
  float **vscAA = NULL;
  int nseq = 0;
  int i, v, idx;

  status = MPI_Unpack (buf, n, pos, &nseq,        1, MPI_INT,   comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  ESL_ALLOC(vscA, sizeof(float) * (M*nseq));
  status = MPI_Unpack (buf, n, pos, vscA, (M*nseq), MPI_FLOAT,  comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  ESL_ALLOC(vscAA, sizeof(float *) * (M));
  idx = 0;
  for(v = 0; v < M; v++) { 
    ESL_ALLOC(vscAA[v], sizeof(float) * nseq);
    for(i = 0; i < nseq; i++)
      vscAA[v][i] = vscA[idx++];
  }
  ESL_DASSERT1((idx == (M*nseq)));

  ESL_ALLOC(vit_cp9scA, sizeof(float) * nseq);
  status = MPI_Unpack (buf, n, pos, vit_cp9scA, nseq, MPI_FLOAT,  comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  ESL_ALLOC(fwd_cp9scA, sizeof(float) * nseq);
  status = MPI_Unpack (buf, n, pos, fwd_cp9scA, nseq, MPI_FLOAT,  comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  ESL_ALLOC(partA,      sizeof(int) * nseq);
  status = MPI_Unpack (buf, n, pos, partA, nseq, MPI_INT,    comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  *ret_vscAA = vscAA;
  *ret_vit_cp9scA = vit_cp9scA;
  *ret_fwd_cp9scA = fwd_cp9scA;
  *ret_partA      = partA;
  *ret_nseq = nseq;
  return eslOK;

 ERROR:
  if(vit_cp9scA != NULL) free(vit_cp9scA);
  if(vit_cp9scA != NULL) free(fwd_cp9scA);
  if(partA      != NULL) free(partA);

  if(vscA  != NULL) free(vscA);
  if(vscAA != NULL) { 
    for(i = 0; i < nseq; i++)
      free(vscAA[i]);
    free(vscAA);
  }
  *ret_vscAA = NULL;
  *ret_nseq = 0;
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
  status = esl_mpi_PackOptSize(comlog->ccom1, -1, MPI_CHAR, comm, &sz); n += sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  status = esl_mpi_PackOptSize(comlog->cdate1,-1, MPI_CHAR, comm, &sz); n += sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  status = esl_mpi_PackOptSize(comlog->ccom2, -1, MPI_CHAR, comm, &sz); n += sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");
  status = esl_mpi_PackOptSize(comlog->cdate2,-1, MPI_CHAR, comm, &sz); n += sz; if (status != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");

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
  status = esl_mpi_PackOpt(comlog->ccom1,  -1, MPI_CHAR, buf, n, position,  comm); if (status != eslOK) ESL_EXCEPTION(eslESYS, "pack failed");
  status = esl_mpi_PackOpt(comlog->cdate1, -1, MPI_CHAR, buf, n, position,  comm); if (status != eslOK) ESL_EXCEPTION(eslESYS, "pack failed");
  status = esl_mpi_PackOpt(comlog->ccom2,  -1, MPI_CHAR, buf, n, position,  comm); if (status != eslOK) ESL_EXCEPTION(eslESYS, "pack failed");
  status = esl_mpi_PackOpt(comlog->cdate2, -1, MPI_CHAR, buf, n, position,  comm); if (status != eslOK) ESL_EXCEPTION(eslESYS, "pack failed");
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

  ESL_ALLOC(comlog, sizeof(ComLog_t));
  status = esl_mpi_UnpackOpt(buf, n, pos, (void**)&(comlog->bcom),   NULL, MPI_CHAR, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = esl_mpi_UnpackOpt(buf, n, pos, (void**)&(comlog->bdate),  NULL, MPI_CHAR, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = esl_mpi_UnpackOpt(buf, n, pos, (void**)&(comlog->ccom1),  NULL, MPI_CHAR, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = esl_mpi_UnpackOpt(buf, n, pos, (void**)&(comlog->cdate1), NULL, MPI_CHAR, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = esl_mpi_UnpackOpt(buf, n, pos, (void**)&(comlog->ccom2),  NULL, MPI_CHAR, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  status = esl_mpi_UnpackOpt(buf, n, pos, (void**)&(comlog->cdate2), NULL, MPI_CHAR, comm); if (status != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  *ret_comlog = comlog;
  return eslOK;

 ERROR:
  if(comlog != NULL) free(comlog);
  *ret_comlog = NULL;
  return status;
}

#endif

