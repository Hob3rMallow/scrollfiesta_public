#ifndef CLOSEUP_WINDOW_INCLUDED
#define CLOSEUP_WINDOW_INCLUDED

/* Fit a close-up crop to the physically matched required pieces, then place an
 * expanded crop without covering stationary strip fragments. JSON plan on stdout. */
int CloseupWindow_run(const char *input);
int CloseupWindow_selftest(void);

#endif
