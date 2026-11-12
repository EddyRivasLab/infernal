/* errors.c
 * 
 * error handling.
 * 
 * Infernal's fatal error messages distinguish between user errors
 * ("failure", with cm_Fail()) and internal faults ("death", with
 * cm_Die()). For now, though, there is no difference between the two
 * functions. Someday we might have cm_Die() print a comforting
 * apology, or provide some help on how to report bugs to us;
 * cm_Fail() might provide some pointers on where to read more
 * documentation.
 * 
 * Based on HMMER3's errors.c:
 * SRE, Fri Jan 12 08:46:02 2007
 */

#include <esl_config.h>
#include <p7_config.h>
#include "config.h"

#include <stdlib.h>
#include <stdio.h>

#include "easel.h"

#include "hmmer.h"

#include "infernal.h"

/* Function:  cm_Die()
 * Incept:    EPN, Fri Jul 27 14:35:44 2007
 *
 * Purpose:   Handle a fatal exception (something that's the system's fault,
 *            including memory allocation failures; or possibly our fault).
 */
void
cm_Die(char *format, ...)
{
  va_list  argp;
                                /* format the error mesg */
  fprintf(stderr, "\nFATAL: ");
  va_start(argp, format);
  vfprintf(stderr, format, argp);
  va_end(argp);
  fprintf(stderr, "\n");
  fflush(stderr);
  exit(1);
}

/* Function:  cm_Fail()
 * Incept:    EPN, Fri Jul 27 14:35:58 2007
 *
 * Purpose:   Handle a user error (something that's the user's fault).
 */
void
cm_Fail(char *format, ...)
{
  va_list  argp;
                                /* format the error mesg */
  fprintf(stderr, "\nError: ");
  va_start(argp, format);
  vfprintf(stderr, format, argp);
  va_end(argp);
  fprintf(stderr, "\n");
  fflush(stderr);
  exit(1);
}


/* Function:  cm_StageVerbose()
 * Incept:    EPN, Sun Oct 11 2026 (w/Claude)
 *
 * Purpose:   Return TRUE if developer diagnostics are on (--devdiag;
 *            see cm_SetStageVerbose()), FALSE if not.
 *
 *            Developer diagnostics -- which alignment engine or band
 *            builder ran, and per-stage timings -- are printed only
 *            when they are on.
 */
static int stage_verbose = FALSE; /* see cm_StageVerbose(), cm_SetStageVerbose() */

int
cm_StageVerbose(void)
{
  return stage_verbose;
}

/* Function:  cm_SetStageVerbose()
 * Incept:    EPN, Tue Oct 20 2026 (w/Claude)
 *
 * Purpose:   Turn developer diagnostics on (<on> = TRUE) or off for the
 *            rest of the run. Programs call this right after option
 *            parsing when --devdiag is used (brief 26_0821-135). The
 *            INFERNAL_STAGE_VERBOSE environment variable that also did
 *            this was removed in brief 26_0821-174.
 */
void
cm_SetStageVerbose(int on)
{
  stage_verbose = on ? TRUE : FALSE;
}
