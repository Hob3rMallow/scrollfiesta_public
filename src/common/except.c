#include "except.h"
#include <stdlib.h>
#include <stdio.h>
#include "ves_omp.h"

const Except_T Arena_Failed = { "Arena allocation failed" };
const Except_T IO_Failed    = { "I/O operation failed" };
const Except_T Timeout      = { "Volume processing timeout" };

EXCEPT_THREAD_LOCAL Except_Frame *Except_stack = NULL;

void Except_raise(const Except_T *e, const char *file, int line)
{
    Except_Frame *p = Except_stack;

    if (p == NULL) {
        fprintf(stderr, "Uncaught exception: %s at %s:%d\n",
                e->reason, file, line);
        abort();
    }
    p->exception = e;
    p->file      = file;
    p->line      = line;
    p->flag      = Except_raised;
    Except_stack = p->prev;
    longjmp(p->env, Except_raised);
}

static int except_nested_control(void)
{
    Except_Frame *initial = Except_stack;
    volatile int failures = 0;
    for (int i = 0; i < 100; i++) {
        volatile int caught = 0;
        TRY {
            TRY { RAISE(Timeout); }
            EXCEPT(IO_Failed) { failures++; }
            END_TRY; /* Reraises into this thread's outer frame. */
        } EXCEPT(Timeout) { caught = 1; }
        END_TRY;
        if (!caught || Except_stack != initial) failures++;
    }
    return failures;
}

int Except_selftest(void)
{
    Except_Frame **slots[4] = {0};
    int count = 1;
    volatile int failures = 0;
    Except_Frame *initial = Except_stack;
    TRY {
        Except_Frame *outer = Except_stack;
#ifdef _OPENMP
#pragma omp parallel num_threads(4)
        {
            int tid = omp_get_thread_num();
            slots[tid] = &Except_stack;
#pragma omp single
            count = omp_get_num_threads();
        }
#else
        slots[0] = &Except_stack;
#endif
        /* Check isolation before throwing: the old global implementation must
         * fail this control without attempting an unsafe cross-thread jump. */
        for (int i = 0; i < count; i++) for (int j = 0; j < i; j++)
            if (slots[i] == slots[j]) failures++;
        if (!failures) {
            int nested_failures = 0, task = 0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(4) reduction(+:nested_failures)
#endif
            for (task = 0; task < 16; task++) nested_failures += except_nested_control();
            failures += nested_failures;
        }
        if (Except_stack != outer) failures++;
    } EXCEPT(Timeout) { failures++; }
    END_TRY;
    if (Except_stack != initial) failures++;
    fprintf(stderr,"  Exception isolation: %d threads, nested rethrow/restoration; %d failures\n",count,(int)failures);
    return failures;
}
