#ifndef CLOSEUP_SEED_INCLUDED
#define CLOSEUP_SEED_INCLUDED

/* Exact, bounded crop of the rendered source+fill+growth mesh. The manifest
 * binds component chart sidecars, the region index, box, margin and output.
 * No fit, welding, coordinate change or chart renumbering is performed. */
int CloseupSeed_write(const char *manifest);
int CloseupSeed_selftest(void);

#endif
