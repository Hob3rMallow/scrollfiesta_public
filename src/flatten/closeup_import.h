#ifndef CLOSEUP_IMPORT_INCLUDED
#define CLOSEUP_IMPORT_INCLUDED

/* Preserve a reviewed reading-canvas mesh, adding only an integer UV translation
 * and renderer companions. No fitting, cuts, welding or face removal. The caller
 * supplies the source's per-vertex chart IDs and a new output stem. */
int CloseupImport_write(const char *mesh, const char *charts, const char *stem, const char *chart_map);
int CloseupImport_selftest(void);

#endif
