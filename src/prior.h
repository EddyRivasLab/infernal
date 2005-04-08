/* prior.h
 * Dirichlet priors for parameterizing a new model.
 * 
 * Original code from Eric Nawrocki; adapted by SRE.
 * SRE, Thu Apr  7 10:16:54 2005
 * SVN $Id$
 */

#define MAXDCHLET     20        /* Maximum number of mixture components */
#define NTRANSSETS    74        /* Number of transition sets (should be 74) */

/* Structure: Prior_t
 * 
 * Dirichlet priors on model parameters. Hacked out of HMMER.
 * This structure contains all the priors - transitions, match bp emissions,
 * singlet emissions, and insert emissions.
 *
 * Fixed-size allocations for convenience; may want to make this
 * dynamic eventually.
 */
struct {
  /* transition priors */
  int    tsetnum;                           /* number of transition sets to read in */
  int    tsetmap[UNIQUESTATES][NODETYPES];  /* tsetmap[a][b] is for transition set from ustate a to node b */
  int    tnalpha[NTRANSSETS];               /* # of parameters in each transition set */
  int    tnum[NTRANSSETS];	            /* # of mixture components for each set */
  double tq[NTRANSSETS][MAXDCHLET];         /* mixture coefficients for each set */
  double t[NTRANSSETS][MAXDCHLET][MAXCONNECT];  /* transition terms for each mixture component */

  /* match base pair priors */
  int    mbpnum;                  /* number of match Dirichlet mixture terms      */
  double mbpq[MAXDCHLET];          /* weights of Dirichlet mixture terms     */
  double mbp[MAXDCHLET][(MAXABET*MAXABET)]; /* match emission terms per mix component */
  
  /* match singlet priors */
  int    mntnum;                  /* number of match Dirichlet mixture terms  */
  double mntq[MAXDCHLET];         /* weights of Dirichlet mixture terms     */
  double mnt[MAXDCHLET][MAXABET]; /* match emission terms per mix component */
  
  /* insert priors */
  int    inum;                   /* number of Dirichlet mixture terms      */
  double iq[MAXDCHLET];          /* weights of Dirichlet mixture terms     */
  double i[MAXDCHLET][MAXABET];  /* insert emission terms per mix component */
} Prior_t;


extern Prior_t *Prior_Read(void);


