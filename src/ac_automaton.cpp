#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <vector>
#include <queue>

using namespace std;

#include "timer.h"
#include "utils.h"
#include "ac_automaton.h"
#include "swsharp/swsharp.h"

// ***************************************************************************
// PUBLIC

extern void automatonCreate(Chain** queries, int querieslen, Seeds* seeds,
    int seedLen, ACNode** automaton);

extern void automatonDelete(ACNode* automaton);

// ***************************************************************************

// ***************************************************************************
// PRIVATE

static void automatonAddWord(ACNode* automaton, char* word, int wordLen, 
    int queryIdx, int location);

static void automatonSetSupply(ACNode* automaton);

static int seedCode(Chain* chain, int pos, int seedLen);

static int seedCode(char* seed, int seedLen);

// ***************************************************************************

// ***************************************************************************
// PUBLIC

extern void automatonCreate(Chain** queries, int queriesLen, Seeds* seeds,
    int seedLen, ACNode** automaton) {

    (*automaton) = new ACNode();

    for (int queryIdx = 0; queryIdx < queriesLen; ++queryIdx) {

        Chain* query = queries[queryIdx];
        int queryLen = chainGetLength(query);

        for (int qstart = 0; qstart < queryLen - seedLen + 1; ++qstart) {

            int code = seedCode(query, qstart, seedLen);
            int size = (*seeds)[code].size();

            for (int i = 0; i < size; ++i) {
                automatonAddWord(*automaton, (*seeds)[code][i], seedLen,
                    queryIdx, qstart);
            }
        }
    }

    automatonSetSupply(*automaton);
}

extern void automatonDelete(ACNode* automaton) {
    queue<ACNode*> nodeQ;

    for (int i = 0; i < 26; ++i) {
        if (automaton->edge[i] && automaton->edge[i] != automaton) {
            nodeQ.push(automaton->edge[i]);
        }
    }

    while (!nodeQ.empty()) {
        ACNode* curr = nodeQ.front();
        // fprintf(stderr, "%p\n", curr);

        nodeQ.pop();

        for(int i = 0; i < 26; ++i) {
            if (curr->edge[i]) {
                nodeQ.push(curr->edge[i]);
            }
        }

        delete curr;
    }

    delete automaton;
}

// ***************************************************************************

// ***************************************************************************
// PRIVATE

static void automatonAddWord(ACNode* automaton, char* word, int wordLen, 
    int queryIdx, int location) {

    ACNode* q = automaton;

    for (int i = 0; i < wordLen; ++i) {
        if (!q->edge[word[i] - 'A']) {
            // create new node
            ACNode* next = new ACNode();
            q->edge[word[i] - 'A'] = next;
        }

        q = q->edge[word[i] - 'A'];
    }

    q->final = 1;

    if (q->positions.size() == 0) {
        q->positions.push_back(seedCode(word, wordLen));
    }

    q->positions.push_back(queryIdx);
    q->positions.push_back(location);
}

static void automatonSetSupply(ACNode* automaton) {
    ACNode* q = automaton;
    automaton->fail = automaton;

    queue<ACNode*> nodeQ;

    for (int i = 0; i < 26; ++i) {
        if (automaton->edge[i]) {
            automaton->edge[i]->fail = automaton;
            nodeQ.push(automaton->edge[i]);
        } else {
            automaton->edge[i] = automaton;
        }
    }

    while (!nodeQ.empty()) {
        q = nodeQ.front();
        nodeQ.pop();

        for (int i = 0; i < 26; ++i) {
            if (!q->edge[i]) {
                continue;
            }

            ACNode* next = q->edge[i];
            nodeQ.push(next);

            ACNode* ft = q->fail;
            while (!ft->edge[i]) {
                ft = ft->fail;
            }

            next->fail = ft->edge[i];

            if (ft->edge[i]->final) {
                next->final = 1;
            }
        }    
    }
}    

static int seedCode(Chain* chain, int pos, int seedLen) {

    int code = 0;
    int start = 5 * (seedLen - 1);

    for (int i = 0; i < seedLen; ++i) {
        code += static_cast<int>(chainGetCode(chain, pos + i))
            << (start - 5 * i);
    }

    return code;
}

static int seedCode(char* seed, int seedLen) {

    int code = 0;
    int start = 5 * (seedLen - 1);

    for (int i = 0; i < seedLen; ++i) {
        code += static_cast<int>(toupper(seed[i] - 'A')) << (start - 5 * i);
    }

    return code;
}

// ***************************************************************************
