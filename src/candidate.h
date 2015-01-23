#ifndef __CANDIDATEH__
#define __CANDIDATEH__

#include "swsharp/swsharp.h"
#include <vector>

#ifdef __cplusplus
extern "C" {
#endif

struct Candidate {
    int score;
    int idx;

    Candidate(int score_, int idx_) : score(score_), idx(idx_) {
    }
};

typedef struct Candidate Candidate;
typedef std::vector<std::vector<Candidate> > Candidates;

extern void candidatesCreate(Candidates** candidates, int len);

extern void candidatesDelete(Candidates* candidates);

#ifdef __cplusplus
}
#endif
#endif  // __CANDIDATEH__