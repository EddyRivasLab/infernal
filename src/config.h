#ifndef CONFIGH_INCLUDED
#define CONFIGH_INCLUDED

/* config.h
 * SRE, Sun Jun  3 20:22:38 2001 [St. Louis]
 * CVS $Id$
 * 
 * Configurable compile-time constants in INFERNAL.
 */

/* RAMLIMIT (in MB) defines how much memory we're
 * allowed to expend on alignment algorithms without
 * switching to more efficient memory forms - e.g.
 * in smallcyk.c
 */
#ifndef RAMLIMIT
#define RAMLIMIT 0
#endif                                           

#endif /* CONFIGH_INCLUDED */

