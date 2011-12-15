/************************************************************
 * @LICENSE@
 ************************************************************/
/* cp9_modelconfig.c
 * EPN, Wed Dec  5 12:56:42 2007
 *
 * Note: all of these functions originated in cp9.c [EPN 02.27.06]
 * 
 * Configuring a CP9 HMM into different global or local modes.
 * Included in this file are functions for configuring HMMs that were
 * built for 'sub CMs'.
 * 
 */

#include "esl_config.h"
#include "p7_config.h"
#include "config.h"

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <math.h>

#include "easel.h"
#include "esl_msa.h"
#include "esl_vectorops.h"

#include "funcs.h"
#include "structs.h"
