/*
 * taucs_stubs.c -- satisfy references in taucs_linsolve to TAUCS subsystems we
 * deliberately do NOT compile (out-of-core factor/solve/IO, AMWB/Vaidya
 * preconditioners). taucs_linsolve references these in branches that only run
 * when the caller passes taucs.ooc=* / taucs.approximate.amwb=* options -- we
 * never do, so these stubs are never executed (except taucs_available_memory_size,
 * which we make report plenty of RAM so the in-core path is chosen).
 *
 * No taucs.h here on purpose: the linker resolves by name, so simplified
 * signatures are fine and avoid any prototype cross-check.
 */

double taucs_available_memory_size(void)
{
    return 8.0e9;   /* report 8 GB -> taucs_linsolve picks the in-core factor */
}

void *taucs_amwb_preconditioner_create(void *A, int a, double b, int c)
{
    (void)A; (void)a; (void)b; (void)c;
    return 0;       /* no AMWB preconditioner */
}

void *taucs_io_create_multifile(char *filename) { (void)filename; return 0; }
void *taucs_io_open_multifile  (char *filename) { (void)filename; return 0; }
int   taucs_io_close (void *f) { (void)f; return -1; }
int   taucs_io_delete(void *f) { (void)f; return -1; }

int taucs_ooc_factor_llt(void *A, void *handle, double memory)
{
    (void)A; (void)handle; (void)memory;
    return -1;      /* TAUCS_ERROR; OOC factor not compiled in */
}

int taucs_ooc_solve_llt(void *L, void *x, void *b)
{
    (void)L; (void)x; (void)b;
    return -1;      /* OOC solve not compiled in */
}
