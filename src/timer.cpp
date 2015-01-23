#include <stdio.h>
#include <stdlib.h>

#include <timer.h>

// ***************************************************************************
// PUBLIC

extern void timerStart(Timeval* timer);

extern long long timerStop(Timeval* timer);

extern void timerPrint(char message[], long long usec);

// ***************************************************************************

// ***************************************************************************
// PRIVATE

// ***************************************************************************



// ***************************************************************************
// PUBLIC

extern void timerStart(Timeval* timer) {
    gettimeofday(timer, NULL);
}

extern long long timerStop(Timeval* timer) {

    Timeval stop;
    gettimeofday(&stop, NULL);

    return ((stop.tv_sec - timer->tv_sec) *
        1000000L + stop.tv_usec) - timer->tv_usec;
}

extern void timerPrint(char message[], long long usec) {
    fprintf(stderr, "[%-25s]: %25lld us\n", message, usec);
}

// ***************************************************************************

// ***************************************************************************
// PRIVATE

// ***************************************************************************
