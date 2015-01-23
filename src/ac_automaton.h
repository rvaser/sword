#ifndef __AC_AUTOMATONHH__
#define __AC_AUTOMATONHH__

#include "utils.h"
#include "swsharp/swsharp.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
    State of the automaton.

    final - denotes if the current state corresponds to kmers
    sup - fail link of AC automaton
    transitions - transition table
    wordLocations - if state is final, then this list stores 
        the locations of the current kmer

*/

struct ACNode {
    int final;
    ACNode* fail;
    ACNode* edge[26];
    std::vector<int> positions;

    ACNode() : final(0), fail(NULL) {
        for (int i = 0; i < 26; ++i) {
            edge[i] = NULL;
        }
    }
};

typedef struct ACNode ACNode;

extern void automatonCreate(Chain** queries, int queriesLen, Seeds* seeds,
    int seedLen, ACNode** automaton);

extern void automatonDelete(ACNode* automaton);

#ifdef __cplusplus    
}
#endif


#endif // __AC_AUTOMATONHH__