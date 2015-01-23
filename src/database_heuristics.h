#ifndef __DATABASE_HASHH__
#define __DATABASE_HASHH__

#include "swsharp/swsharp.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void* databaseIndicesCreate(char* databasePath, char* queryPath, int seedLen,
    int maxCandidates, int hitThreshold, Scorer* scorer, int threadLen);

extern void databaseIndicesDelete(void* indices_);

extern int* filteredDatabaseCreate(Chain*** filteredDatabase,
    int* filteredDatabaseLen, void* indices_, int queryIdx,
    Chain** database, int databaseLen, int returnUsed);

extern void filteredDatabaseDelete(Chain** filteredDatabase);

#ifdef __cplusplus
}
#endif
#endif  // __DATABASE_HASHH__
