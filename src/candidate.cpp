#include <stdio.h>
#include <stdlib.h>
#include <vector>

using namespace std;

#include "swsharp/swsharp.h"
#include "candidate.h"

// ***************************************************************************
// PUBLIC

extern void candidatesCreate(Candidates** candidates, int len);

extern void candidatesDelete(Candidates* candidates);

// ***************************************************************************

// ***************************************************************************
// PRIVATE

// ***************************************************************************



// ***************************************************************************
// PUBLIC

extern void candidatesCreate(Candidates** candidates, int len) {
    vector<Candidate> vc;
    (*candidates) = new Candidates(len, vc);
}

extern void candidatesDelete(Candidates* candidates) {
    for (unsigned int i = 0; i < candidates->size(); ++i) {
        vector<Candidate>().swap((*candidates)[i]);
    }
    candidates->clear();
    delete candidates;
}

// ***************************************************************************

// ***************************************************************************
// PRIVATE

// ***************************************************************************