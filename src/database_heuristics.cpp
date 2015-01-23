#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <vector>
#include <fstream>
#include <iostream>
#include <utility>
#include <algorithm>
#include <functional>

using namespace std;

#include "timer.h"
#include "utils.h"
#include "candidate.h"
#include "ac_automaton.h"
#include "database_heuristics.h"
#include "swsharp/swsharp.h"

#define ASSERT(expr, fmt, ...)\
    do {\
        if (!(expr)) {\
            fprintf(stderr, "[ERROR]: " fmt "\n", ##__VA_ARGS__);\
            exit(-1);\
        }\
    } while (0)

struct Sequence {
    int idx;
    int len;

    Sequence(int idx_, int len_) :
        idx(idx_), len(len_) {
    }
};

static bool sequenceByLen(const Sequence& left, const Sequence& right) {
    return left.len < right.len;
}

struct ThreadData{
    int threadIdx;
    Chain** queries;
    vector<Sequence>* qsequences;
    Chain** database;
    vector<Sequence>* tsequencesShort;
    vector<int>* tsegmentsShort;
    vector<Sequence>* tsequencesLong;
    vector<int>* tsegmentsLong;
    Seeds* seeds;
    int seedLen;
    Candidates* candidates;
    Mutex* candidatesMutex;
    int maxCandidates;

    ThreadData(int threadIdx_, Chain** queries_, vector<Sequence>* qsequences_,
            Chain** database_, vector<Sequence>* tsequencesShort_,vector<int>* tsegmentsShort_,
            vector<Sequence>* tsequencesLong_, vector<int>* tsegmentsLong_, Seeds* seeds_,
            int seedLen_, Candidates* candidates_, Mutex* candidatesMutex_,
            int maxCandidates_) :
        threadIdx(threadIdx_), queries(queries_), qsequences(qsequences_), database(database_),
        tsequencesShort(tsequencesShort_), tsegmentsShort(tsegmentsShort_),
        tsequencesLong(tsequencesLong_), tsegmentsLong(tsegmentsLong_), seeds(seeds_),
        seedLen(seedLen_), candidates(candidates_), candidatesMutex(candidatesMutex_),
        maxCandidates(maxCandidates_) {
    }
};

typedef struct ThreadData ThreadData;

static bool candidateByScore(const Candidate& left, const Candidate& right) {
    return left.score > right.score;
}

// ***************************************************************************
// PUBLIC

extern void* databaseIndicesCreate(char* databasePath, char* queryPath, int seedLen,
    int maxCandidates, int hitThreshold, Scorer* scorer, int threadLen);

extern void databaseIndicesDelete(void* indices_);

extern int* filteredDatabaseCreate(void* indices_, int queryIdx,
    Chain** database, int databaseLen, Chain*** filteredDatabase,
    int* filteredDatabaseLen, int returnUsed);

extern void filteredDatabaseDelete(Chain** filteredDatabase);

// ***************************************************************************

// ***************************************************************************
// PRIVATE

static void preprocQueries(vector<Sequence>& qsequences, Chain** queries,
    int queriesLen);

static void preprocDatabase(vector<Sequence>& tsequencesShort, vector<int>& tsegmentsShort,
    vector<Sequence>& tsequencesLong, vector<int>& tsegmentsLong, Chain** database,
    int databaseLen, int databaseStart, int threadLen);

static void scoreSequences(int threadIdx, Chain** queries, vector<Sequence>* qsequences,
    Chain** database, vector<Sequence>* tsequences, vector<int>* tsegments, Seeds* seeds,
    int seedLen, Candidates* candidates, Mutex* candidatesMutex, int maxCandidates);

static void* findCandidates(void* params);

// ***************************************************************************



// ***************************************************************************
// PUBLIC

extern void* databaseIndicesCreate(char* databasePath, char* queryPath, int seedLen,
    int maxCandidates, int hitThreshold, Scorer* scorer, int threadLen) {

    Timeval timer;
    timerStart(&timer);

    Chain** queries = NULL;
    int queriesLen = 0;
    readFastaChains(&queries, &queriesLen, queryPath);

    vector<Sequence> qsequences;
    preprocQueries(qsequences, queries, queriesLen);

    Seeds* seeds = NULL;
    if (seedLen == 3) {
        seedsCreateLong(&seeds, seedLen, hitThreshold, scorer);
    } else {
        // only 1 substitution is allowed
        seedsCreateShort(&seeds, seedLen, hitThreshold, scorer);
    }

    Candidates* candidates = NULL;
    candidatesCreate(&candidates, queriesLen);

    Data* indices = NULL;
    dataCreate(&indices, queriesLen);

    Chain** database = NULL;
    int databaseLen = 0;
    int databaseStart = 0;

    FILE* handle;
    int serialized;

    readFastaChainsPartInit(&database, &databaseLen, &handle, &serialized, databasePath);

    Mutex candidatesMutex;
    mutexCreate(&candidatesMutex);

    while (1) {

        int status = 1;

        status &= readFastaChainsPart(&database, &databaseLen, handle,
            serialized, 1000000000); // ~1GB

        vector<Sequence> tsequencesShort;
        vector<int> tsegmentsShort;
        vector<Sequence> tsequencesLong;
        vector<int> tsegmentsLong;
        preprocDatabase(tsequencesShort, tsegmentsShort, tsequencesLong, tsegmentsLong,
            database, databaseLen, databaseStart, threadLen);

        ThreadPoolTask** threadTasks = new ThreadPoolTask*[threadLen];

        for (int i = 0; i < threadLen; ++i) {

            ThreadData* threadData = new ThreadData(i, queries, &qsequences, database,
                &tsequencesShort, &tsegmentsShort, &tsequencesLong, &tsegmentsLong,
                seeds, seedLen, candidates, &candidatesMutex, maxCandidates);

            threadTasks[i] = threadPoolSubmit(findCandidates, static_cast<void*>(threadData));
        }

        for (int i = 0; i < threadLen; ++i) {
            threadPoolTaskWait(threadTasks[i]);
            threadPoolTaskDelete(threadTasks[i]);
        }

        delete[] threadTasks;

        if (status == 0) {
            break;
        }

        for (int i = databaseStart; i < databaseLen; ++i) {
            chainDelete(database[i]);
            database[i] = NULL;
        }

        databaseStart = databaseLen;
    }

    for (int i = 0; i < queriesLen; ++i) {
        (*indices)[i].reserve((*candidates)[i].size());

        for (int j = 0; j < (int) (*candidates)[i].size(); ++j) {
            (*indices)[i].push_back((*candidates)[i][j].idx);
        }

        sort((*indices)[i].begin(), (*indices)[i].end());
    }

    mutexDelete(&candidatesMutex);

    fclose(handle);

    deleteFastaChains(database, databaseLen);
    deleteFastaChains(queries, queriesLen);

    candidatesDelete(candidates);

    seedsDelete(seeds);

    timerPrint("HeuristicsTotal", timerStop(&timer));

    return static_cast<void*>(indices);
}

extern void databaseIndicesDelete(void* indices_) {
    Data* indices = static_cast<Data*>(indices_);
    dataDelete(indices);
}

extern int* filteredDatabaseCreate(Chain*** filteredDatabase,
    int* filteredDatabaseLen, void* indices_, int queryIdx,
    Chain** database, int databaseLen, int returnUsed) {

    Data* indices = static_cast<Data*>(indices_);

    int* usedIndices = NULL;
    int usedIndicesLen = 0;

    int databaseEnd = databaseLen - 1;
    unsigned int j;

    for (j = 0; j < (*indices)[queryIdx].size(); ++j) {
        if ((*indices)[queryIdx][j] > databaseEnd) break;
    }

    usedIndicesLen = j;

    if (usedIndicesLen == 0) {
        *filteredDatabase = NULL;
        *filteredDatabaseLen = 0;
    } else {
        *filteredDatabase = new Chain*[usedIndicesLen];
        *filteredDatabaseLen = usedIndicesLen;

        for (int i = 0; i < usedIndicesLen; ++i) {
            (*filteredDatabase)[i] = database[(*indices)[queryIdx][i]];
        }

        if (returnUsed) {
            usedIndices = static_cast<int*>(malloc(usedIndicesLen * sizeof(*usedIndices)));

            for (int i = 0; i < usedIndicesLen; ++i) {
                usedIndices[i] = (*indices)[queryIdx][i];
            }
        }

        vector<int> temp(
            (*indices)[queryIdx].begin() + usedIndicesLen,
            (*indices)[queryIdx].end());

        (*indices)[queryIdx].swap(temp);
    }

    return usedIndices;
}

extern void filteredDatabaseDelete(Chain** filteredDatabase) {
    delete[] filteredDatabase;
}

// ***************************************************************************

// ***************************************************************************
// PRIVATE

static void preprocQueries(vector<Sequence>& qsequences, Chain** queries,
    int queriesLen) {

    int totalLen = 0;

    qsequences.reserve(queriesLen);

    for (int i = 0; i < queriesLen; ++i) {
        int len = chainGetLength(queries[i]);

        qsequences.emplace_back(i, len);
        totalLen += len;
    }

    sort(qsequences.begin(), qsequences.end(), sequenceByLen);
}

static void preprocDatabase(vector<Sequence>& tsequencesShort, vector<int>& tsegmentsShort,
    vector<Sequence>& tsequencesLong, vector<int>& tsegmentsLong, Chain** database,
    int databaseLen, int databaseStart, int threadLen) {

    int totalShortLen = 0;
    int totalLongLen = 0;

    for (int i = databaseStart; i < databaseLen; ++i) {
        int len = chainGetLength(database[i]);

        if (len > 2000) {
            tsequencesLong.emplace_back(i, len);
            totalLongLen += len;
        } else {
            tsequencesShort.emplace_back(i, len);
            totalShortLen += len;
        }
    }

    // distribute short sequences
    sort(tsequencesShort.begin(), tsequencesShort.end(), sequenceByLen);

    tsegmentsShort.push_back(0);

    int segmentMaxLen = totalShortLen / threadLen;
    int segmentLen = 0;

    for (int i = 0; i < (int) tsequencesShort.size(); ++i) {

        segmentLen += tsequencesShort[i].len;

        if (segmentLen > segmentMaxLen) {
            tsegmentsShort.push_back(i);
            segmentLen = 0;

            if ((int) tsegmentsShort.size() == threadLen) {
                break;
            }
        }
    }

    while ((int) tsegmentsShort.size() != threadLen) {
        tsegmentsShort.push_back(tsegmentsShort.back());
    }

    tsegmentsShort.push_back((int) tsequencesShort.size());

    fprintf(stderr, "[Short]\n");
    for (int i = 0; i < (int) tsegmentsShort.size(); ++i) {
        fprintf(stderr, "[%d] %d\n", i, tsegmentsShort[i]);
    }

    // distribute long sequences
    sort(tsequencesLong.begin(), tsequencesLong.end(), sequenceByLen);

    tsegmentsLong.push_back(0);

    if (threadLen > 1) {
        int workingLen = threadLen;
        int stopSize = threadLen;

        if (threadLen > 3) {
            workingLen -= 2;
            stopSize -= 1;
            tsegmentsLong.push_back(0);  
        }

        segmentMaxLen = totalLongLen / (workingLen);
        segmentLen = 0;

        for (int i = 0; i < (int) tsequencesLong.size(); ++i) {

            segmentLen += tsequencesLong[i].len;

            if (segmentLen > segmentMaxLen) {
                tsegmentsLong.push_back(i);
                segmentLen = 0;

                if ((int) tsegmentsLong.size() == stopSize) {
                    break;
                }
            }
        }

        while ((int) tsegmentsLong.size() != stopSize) {
            tsegmentsLong.push_back(tsegmentsLong.back());
        }

        if (threadLen > 3) {
            tsegmentsLong.push_back((int) tsequencesLong.size());
        }
    }

    tsegmentsLong.push_back((int) tsequencesLong.size());

    fprintf(stderr, "[Long]\n");
    for (int i = 0; i < (int) tsegmentsLong.size(); ++i) {
        fprintf(stderr, "[%d] %d\n", i, tsegmentsLong[i]);
    }
}

static void scoreSequences(int threadIdx, Chain** queries, vector<Sequence>* qsequences,
    Chain** database, vector<Sequence>* tsequences, vector<int>* tsegments, Seeds* seeds,
    int seedLen, Candidates* candidates, Mutex* candidatesMutex, int maxCandidates) {

    int targetStart = tsegments->at(threadIdx);
    int targetEnd = tsegments->at(threadIdx + 1);

    if (targetEnd - targetStart == 0) return;

    int queriesLen = (*qsequences).size();
    int maxTargetLen = (*tsequences)[targetEnd - 1].len;
    int groups = 0;

    Candidates* candidatesPart = NULL;
    candidatesCreate(&candidatesPart, queriesLen);

    ACNode* automaton = NULL;

    for (int queryIdx = 0; queryIdx < queriesLen;) {

        ++groups;

        // create automaton
        int scoresLen = 0;
        int groupLen = 0;

        for (int i = queryIdx; i < queriesLen; ++i) {

            int len = (*qsequences)[i].len + maxTargetLen -
                2 * seedLen + 1;

            if (scoresLen + len > 250000) { // ~0.5MB
                break;
            }

            scoresLen += len;
            ++groupLen;
        }

        Chain** queriesPart = new Chain*[groupLen];

        for (int i = 0; i < groupLen; ++i) {
            queriesPart[i] = queries[(*qsequences)[queryIdx + i].idx];
        }

        automatonCreate(queriesPart, groupLen, seeds, seedLen, &automaton);

        vector<int> dlens(groupLen, 0);
        vector<int> dstarts(groupLen + 1, 0);

        vector<unsigned short> max(groupLen, 0);
        vector<unsigned short> min(groupLen, 0);

        vector<int> found(groupLen, 0);

        vector<unsigned short> scores(scoresLen, 0);

        mutexLock(candidatesMutex);

        for (int i = 0; i < groupLen; ++i) {
            int qidx = (*qsequences)[queryIdx + i].idx;

            found[i] = (*candidates)[qidx].size();

            min[i] = found[i] > 0 ? (*candidates)[qidx].back().score : 100000000;
        }

        mutexUnlock(candidatesMutex);

        for (int targetIdx = targetStart; targetIdx < targetEnd; ++targetIdx) {

            ACNode* state = automaton;

            Chain* target = database[(*tsequences)[targetIdx].idx];
            int targetLen = chainGetLength(target);
            const char* tcodes = chainGetCodes(target);

            for (int i = 0; i < groupLen; ++i) {
                dlens[i] = (*qsequences)[queryIdx + i].len + targetLen - 2 * seedLen + 1;
                dstarts[i + 1] = dstarts[i] + dlens[i];
            }

            // find hits
            for (int k = 0; k < targetLen; ++k) {
                int c = tcodes[k];

                while (!state->edge[c]) {
                    state = state->fail;
                }
                if (state->edge[c] == state) continue;

                state = state->edge[c];

                if (state->final) {
                    int tstart = k - seedLen + 1;

                    for (int i = 1; i < (int) state->positions.size(); i += 2) {
                        int idx = state->positions[i];
                        int d = (tstart - state->positions[i + 1] + dlens[idx]) %
                            dlens[idx] + dstarts[idx];

                        ++scores[d];

                        if (max[idx] < scores[d]) {
                            max[idx] = scores[d];
                        }
                    }
                }
            }

            // create canidates
            for (int i = 0; i < (int) groupLen; ++i) {
                if (max[i] == 0) {
                    continue;
                }

                int qidx = (*qsequences)[queryIdx + i].idx;
                int flag = (int) (*candidatesPart)[qidx].size() < maxCandidates &&
                    found[i] < maxCandidates;

                if (flag || max[i] >= min[i]) {
                    (*candidatesPart)[qidx].emplace_back(max[i], (*tsequences)[targetIdx].idx);

                    if (min[i] > max[i]) {
                        min[i] = max[i];
                    }
                }
            }

            // clear hits, reset max and scores
            for (int i = 0; i < (int) groupLen; ++i) {
                if (max[i] == 0) {
                    continue;
                }

                fill_n(scores.begin() + dstarts[i], dlens[i], 0);
            }

            fill(max.begin(), max.end(), 0);
        }

        // sort and pick top candidates
        mutexLock(candidatesMutex);

        for (int i = 0; i < (int) groupLen; ++i) {
            int qidx = (*qsequences)[queryIdx + i].idx;

            (*candidates)[qidx].insert((*candidates)[qidx].end(),
                (*candidatesPart)[qidx].begin(),
                (*candidatesPart)[qidx].end());

            stable_sort(
                (*candidates)[qidx].begin(),
                (*candidates)[qidx].end(),
                candidateByScore);

            if ((int) (*candidates)[qidx].size() > maxCandidates) {
                vector<Candidate> temp(
                    (*candidates)[qidx].begin(),
                    (*candidates)[qidx].begin() + maxCandidates);

                (*candidates)[qidx].swap(temp);
            }

            vector<Candidate>().swap((*candidatesPart)[qidx]);
        }

        mutexUnlock(candidatesMutex);

        // delete automaton
        automatonDelete(automaton);

        delete[] queriesPart;

        queryIdx += groupLen;
    }

    candidatesDelete(candidatesPart);
}

static void* findCandidates(void* params) {

    Timeval threadTimer, shortTimer, longTimer;
    timerStart(&threadTimer);

    ThreadData* threadData = static_cast<ThreadData*>(params);

    // score short
    timerStart(&shortTimer);

    scoreSequences(threadData->threadIdx, threadData->queries, threadData->qsequences,
        threadData->database, threadData->tsequencesShort, threadData->tsegmentsShort,
        threadData->seeds, threadData->seedLen, threadData->candidates,
        threadData->candidatesMutex, threadData->maxCandidates);

    long long shortTotal = timerStop(&shortTimer);

    // score long
    timerStart(&longTimer);

    scoreSequences(threadData->threadIdx, threadData->queries, threadData->qsequences,
        threadData->database, threadData->tsequencesLong, threadData->tsegmentsLong,
        threadData->seeds, threadData->seedLen, threadData->candidates,
        threadData->candidatesMutex, threadData->maxCandidates);

    long long longTotal = timerStop(&longTimer);

    fprintf(stderr, "ThreadIdx: [%d]\n", threadData->threadIdx);
    timerPrint("ThreadTotal", timerStop(&threadTimer));
    timerPrint("LongSolver", longTotal);
    timerPrint("ShortSolver", shortTotal);

    delete threadData;

    return NULL;
}

// ***************************************************************************
