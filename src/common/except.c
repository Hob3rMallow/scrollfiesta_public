#include "except.h"
#include <stdlib.h>
#include <stdio.h>

const Except_T Arena_Failed = { "Arena allocation failed" };
const Except_T IO_Failed    = { "I/O operation failed" };
const Except_T Timeout      = { "Volume processing timeout" };

Except_Frame *Except_stack = NULL;

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
