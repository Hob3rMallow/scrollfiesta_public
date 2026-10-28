#ifndef CLOSEUP_CTTRACK_INCLUDED
#define CLOSEUP_CTTRACK_INCLUDED

/* Experimental native CT tracking audit. It reads observed surface vertices
 * at a fixed world-z plane, then follows CT alone toward decreasing z. The
 * withheld geometry is used only for reporting, never during tracking.
 * No surface is admitted or modified by this audit. box is z0,z1,y0,y1,x0,x1. */
int CloseupCtTrack_audit(const char *surface_stem, const char *raw,
                         const char *out_stem, const double box[6], double seed_z);
int CloseupCtTrack_selftest(void);

#endif
