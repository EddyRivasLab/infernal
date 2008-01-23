/* scancyk.c
 * SRE, Thu May  2 11:50:48 2002 [AA 3050 SFO->STL]
 * SVN $Id$
 * 
 * CYK alignment: multihit, local, database scanning mode.
 * [xref STL6 p47]
 * 
 ***************************************************************** 
 * @LICENSE@
 ***************************************************************** 
 */

#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>

#include "easel.h"
#include "esl_sqio.h"
#include "esl_vectorops.h"

#include "funcs.h"
#include "structs.h"
