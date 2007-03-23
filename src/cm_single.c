/************************************************************
 * @LICENSE@
 ************************************************************/

/* cm_single.c
 * EPN, Wed Mar 21 17:25:21 2007
 * 
 * Functions to support building single sequence CMs
 * in cmbuild. 
 * 
 */

#include "config.h"
#include "squidconf.h"

#include <stdio.h>
#include <math.h>
#include <float.h>

#include "easel.h"		
#include "squid.h"
#include "structs.h"
#include "funcs.h"

/* Function: DivideMSA2SingleMSAs()
 * EPN, Wed Mar 21 17:26:39 2007
 * 
 * Purpose:  Given an MSA, divide it into multiple MSAs, each with
 *           a single sequence, each to be used to build a CM from.
 *           If do_sall is TRUE, number of new MSAs is number of
 *           seqs in input MSA. If do_sall is FALSE, we choose 
 *           nrep representative seqs from the MSA and return 
 *           nrep new MSAs.
 * 
 * Args:    
 * MSA           *mmsa - the master MSA
 * int         do_sall - TRUE to build MSA from each seq in mMSA 
 * int            nrep - number of representative seqs to build MSAs for
 * int       *ret_nmsa - number of MSAs in ret_MSA
 * MSA     ***ret_smsa - new single sequence MSAs
 *           
 * Return:   ret_MSA (alloc'ed here) and ret_nMSA
 */
int 
DivideMSA2SingleMSAs(MSA *mmsa, int do_sall, int nrep, int *ret_nmsa, MSA ***ret_smsa)
{
  int *useme;  /* [0..mmsa->nseq-1], 1 to build single seq MSA from seq i, 0 not to */
  int nmsa;    /* number of single seq MSAs to create */
  MSA **smsa;  /* the single seq MSAs */
  int i;       /* counter over sequences */
  int m;       /* counter over single MSAs */
  char buffer[50];
  int n;
  nmsa = 0;
  useme = MallocOrDie(sizeof(int) * mmsa->nseq);
  if(do_sall)
    for(i = 0; i < mmsa->nseq; i++)
      {
	useme[i] = TRUE;
	nmsa++;
      }

  smsa = MallocOrDie(sizeof(MSA *) * nmsa);
  m = 0;
  for(i = 0; i < mmsa->nseq; i++)
    if(useme[i])
      {
	smsa[m]             = MSAAlloc(1, mmsa->alen);
	smsa[m]->aseq[0]    = sre_strdup(mmsa->aseq[i], -1);
	if(mmsa->desc != NULL) smsa[m]->desc = sre_strdup(mmsa->desc, -1);
	smsa[m]->ss_cons    = sre_strdup(mmsa->ss_cons, -1);
	smsa[m]->nseq       = 1;
	n = sprintf (buffer, ".%d", (m+1));
	smsa[m]->name       = sre_strdup(mmsa->name, -1);
	sre_strcat(&smsa[m]->name, -1, buffer, (n+1));
	printf("smsa[m:%d]->name: %s\n", m, smsa[m]->name);
	m++;
      }

  *ret_nmsa = nmsa;
  *ret_smsa = smsa;

  return eslOK;
}
