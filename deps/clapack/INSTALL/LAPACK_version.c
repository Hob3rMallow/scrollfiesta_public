#include "blaswrap.h"
#include <stdio.h>
#include "f2c.h"

/* Main program */ int MAIN__(void)
{
    /* Local variables */
    static integer patch, major, minor;
    extern /* Subroutine */ int ilaver_(integer *, integer *, integer *);
    ilaver_(&major, &minor, &patch);
    printf("LAPACK %d.%d.%d \n",major,minor,patch);
    return 0;
} /* MAIN__   
   Main program alias */ int lapack_version__ () { MAIN__ (); return 0; }


