/* dsecnd_NONE.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

doublereal dsecnd_(void)
{
    /* System generated locals */
    doublereal ret_val;


/*  -- LAPACK auxiliary routine (version 3.1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     February 2007 */

/*  Purpose */
/*  ======= */

/*  DSECND returns nothing instead of returning the user time for a process in seconds. */
/*  If you are using that routine, it means that neither EXTERNAL ETIME, */
/*  EXTERNAL ETIME_, INTERNAL ETIME, INTERNAL CPU_TIME is available  on */
/*  your machine. */

/* ===================================================================== */

    ret_val = 0.;
    return ret_val;

/*     End of DSECND */

} /* dsecnd_ */
