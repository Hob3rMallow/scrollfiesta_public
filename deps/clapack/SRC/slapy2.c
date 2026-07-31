#include "blaswrap.h"
#include "f2c.h"

doublereal slapy2_(real *x, real *y)
{
/*  -- LAPACK auxiliary routine (version 3.1) --   
       Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
       November 2006   


    Purpose   
    =======   

    SLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary   
    overflow.   

    Arguments   
    =========   

    X       (input) REAL   
    Y       (input) REAL   
            X and Y specify the values x and y.   

    ===================================================================== */
    /* System generated locals */
    real ret_val, r__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    static real w, z__, xabs, yabs;



    xabs = dabs(*x);
    yabs = dabs(*y);
    w = dmax(xabs,yabs);
    z__ = dmin(xabs,yabs);
    if (z__ == 0.f) {
	ret_val = w;
    } else {
/* Computing 2nd power */
	r__1 = z__ / w;
	ret_val = w * sqrt(r__1 * r__1 + 1.f);
    }
    return ret_val;

/*     End of SLAPY2 */

} /* slapy2_ */
