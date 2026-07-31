/* second_EXT_ETIME.f -- translated by f2c (version 20061008).
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

doublereal second_(void)
{
    /* System generated locals */
    real ret_val;

    /* Local variables */
    static real t1;
    extern doublereal etime_(real *);
    static real tarray[2];


/*  -- LAPACK auxiliary routine (version 3.1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     February 2007 */

/*  Purpose */
/*  ======= */

/*  SECOND returns the user time for a process in seconds. */
/*  This version gets the time from the EXTERNAL system function ETIME. */

/* ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    t1 = etime_(tarray);
    ret_val = tarray[0];
    return ret_val;

/*     End of SECOND */

} /* second_ */
