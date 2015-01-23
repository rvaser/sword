#ifndef __UTILSH__
#define __UTILSH__

#include "swsharp/swsharp.h"
#include <vector>

#ifdef __cplusplus
extern "C" {
#endif

typedef std::vector<std::vector<int> > Data;
typedef std::vector<std::vector<char*> > Seeds;

extern void seedsCreateLong(Seeds** seeds, int seedLen, int hitThreshold, Scorer* scorer);

extern void seedsCreateShort(Seeds** seeds, int seedLen, int hitThreshold, Scorer* scorer);

extern void seedsDelete(Seeds* seeds);

extern void seedScoresCreate(int** seedScores, int* seedScoresLen, int seedLen,
    Scorer* scorer);

extern void seedScoresDelete(int* seedScores);

extern void dataCreate(Data** data, int len);

extern void dataDelete(Data* data);


#ifdef __cplusplus
}
#endif
#endif  // __UTILSH__
