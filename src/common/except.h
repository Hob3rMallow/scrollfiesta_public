#ifndef EXCEPT_INCLUDED
#define EXCEPT_INCLUDED

#include <setjmp.h>
#include <stdio.h>

/* Exception type - Hanson Ch.4 */
typedef struct {
    const char *reason;
} Except_T;

/* Global exception instances */
extern const Except_T Arena_Failed;
extern const Except_T IO_Failed;
extern const Except_T Timeout;

/* Exception frame - linked stack */
typedef struct Except_Frame {
    struct Except_Frame *prev;
    jmp_buf             env;
    const Except_T     *exception;
    const char         *file;
    int                 line;
    int                 flag;
} Except_Frame;

enum {
    Except_entered   = 0,
    Except_raised    = 1,
    Except_handled   = 2,
    Except_finalized = 3
};

extern Except_Frame *Except_stack;

void Except_raise(const Except_T *e, const char *file, int line);

#define RAISE(e) Except_raise(&(e), __FILE__, __LINE__)

/* Hanson Ch.4 macros - careful brace matching */
#define TRY do { \
    volatile int Except_flag; \
    Except_Frame Except_frame; \
    Except_frame.prev = Except_stack; \
    Except_stack = &Except_frame; \
    Except_flag = setjmp(Except_frame.env); \
    if (Except_flag == Except_entered) {

#define EXCEPT(e) \
        if (Except_flag == Except_entered) \
            Except_stack = Except_stack->prev; \
    } else if (Except_frame.exception == &(e)) { \
        Except_flag = Except_handled;

#define FINALLY \
        if (Except_flag == Except_entered) \
            Except_stack = Except_stack->prev; \
    } { \
        if (Except_flag == Except_entered) \
            Except_flag = Except_finalized;

#define END_TRY \
        if (Except_flag == Except_entered) \
            Except_stack = Except_stack->prev; \
    } \
    if (Except_flag == Except_raised) \
        Except_raise(Except_frame.exception, \
                     Except_frame.file, Except_frame.line); \
} while (0)

#define RETURN \
    do { \
        Except_stack = Except_stack->prev; \
        return; \
    } while (0)

#endif
