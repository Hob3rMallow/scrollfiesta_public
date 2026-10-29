#ifndef ASM_READING_ORDER_INCLUDED
#define ASM_READING_ORDER_INCLUDED
#include "asm_types.h"
#include "asm_axis.h"

/* Read-only physical navigation. Does not change placements, connect sheets,
 * paint missing material or certify CT identity. Relative turn coordinates
 * are meaningful only inside their reported, consistent evidence group. */
int AsmReadingOrder_write(const AsmRun *run,const AsmAxis *axis,const char *dir);
/* Approximate presentation coordinates, NOT certified cross-piece joins.
 * Source phase is lifted only over measured passing seams; the free gauge
 * of each piece is set by its area-weighted radius / global pitch. */
typedef struct AsmRibbonKey {
    double turn_lo,turn_hi,angular_gradient[2];
    int source_supported;
} AsmRibbonKey;
int AsmReadingOrder_ribbon_keys(const AsmRun *run,const AsmAxis *axis,AsmRibbonKey *keys);
int AsmReadingOrder_selftest(void);
#endif
