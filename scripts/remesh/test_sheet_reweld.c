/*
 * test_sheet_reweld.c -- unit-test driver for the phase-2 sheet-correspondence
 * permissive re-weld (src/remesh/sheet_reweld.c). Runs SheetReweld_selftest:
 *   [A] connected-component sheet labeling counts,
 *   [B1] one clear cross-seam sheet pair is confirmed,
 *   [B2] an ambiguous sheet (overlap split two ways) is rejected,
 *   [C] a confirmed pair actually welds on a restricted (two-sheet) cloud.
 */
#include <stdio.h>
#include "../../src/remesh/sheet_reweld.h"

int main(void)
{
    printf("== sheet_reweld self-test ==\n");
    return SheetReweld_selftest();
}
