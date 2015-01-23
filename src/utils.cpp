#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include <vector>
#include <cstring>

#include "utils.h"
#include "swsharp/swsharp.h"

using namespace std;

#define SEED_IDX_LEN(n) ((n) == 3 ? 26426 : ((n) == 4 ? 845626 : 27060026))

#define AA 20

static const char AMINO_ACIDS[] = {
    'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
    'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '\0'
};

// ***************************************************************************
// PUBLIC

extern void seedsCreateLong(Seeds** seeds, int seedLen, int hitThreshold, Scorer* scorer);

extern void seedsCreateShort(Seeds** seeds, int seedLen, int hitThreshold, Scorer* scorer);

extern void seedsDelete(Seeds* seeds);

extern void seedScoresCreate(int** seedScores, int* seedScoresLen, int seedLen,
    Scorer* scorer);

extern void seedScoresDelete(int* seedScores);

extern void dataCreate(Data** data, int len);

extern void dataDelete(Data* data);

// ***************************************************************************

// ***************************************************************************
// PRIVATE

static int seedCode(char* seed, int seedLen);

static int seedCode(vector<char>* seed);

static void kmersCreate(char*** kmers, int* kmersLen, int seedLen);

static int kmersCreateRec(char*** kmers, int kmerIdx, vector<char>* kmer,
    int n, int seedLen);

static void kmersDelete(char** kmers, int kmersLen);

static int permutation(char* kmer1, char* kmer2, int seedLen, Scorer* scorer);

static int seedScore(vector<char>* seed, Scorer* scorer);

static void seedScoresCreateRec(int** seedScores, vector<char>* seed, int n, Scorer* scorer);

// ***************************************************************************



// ***************************************************************************
// PUBLIC

extern void seedsCreateLong(Seeds** seeds, int seedLen, int hitThreshold, Scorer* scorer) {

    vector<char*> vc;
    (*seeds) = new Seeds(SEED_IDX_LEN(seedLen), vc);

    char** kmers = NULL;
    int kmersLen = 0;
    kmersCreate(&kmers, &kmersLen, seedLen);

    int idxi, idxj;
    int codei, codej;

    for (int i = 0; i < kmersLen; ++i) {
        codei = seedCode(kmers[i], seedLen);
        idxi = (**seeds)[codei].size();

        (**seeds)[codei].push_back(new char[seedLen + 1]);
        strcpy((**seeds)[codei][idxi], kmers[i]);
    }

    if (hitThreshold != -1) {
        for (int i = 0; i < kmersLen; ++i) {
            for (int j = i; j < kmersLen; ++j) {
                if (i == j)  continue;

                if (permutation(kmers[i], kmers[j], seedLen, scorer) >= hitThreshold) {
                    codei = seedCode(kmers[i], seedLen);
                    idxi = (**seeds)[codei].size();

                    (**seeds)[codei].push_back(new char[seedLen + 1]);
                    strcpy((**seeds)[codei][idxi], kmers[j]);

                    codej = seedCode(kmers[j], seedLen);
                    idxj = (**seeds)[codej].size();

                    (**seeds)[codej].push_back(new char[seedLen + 1]);
                    strcpy((**seeds)[codej][idxj], kmers[i]);
                }
            }
        }
    }

    kmersDelete(kmers, kmersLen);
}

extern void seedsCreateShort(Seeds** seeds, int seedLen, int hitThreshold, Scorer* scorer) {

    vector<char*> vc;
    (*seeds) = new Seeds(SEED_IDX_LEN(seedLen), vc);

    char** kmers = NULL;
    int kmersLen = 0;
    kmersCreate(&kmers, &kmersLen, seedLen);

    int idxi, codei;

    for (int i = 0; i < kmersLen; ++i) {
        codei = seedCode(kmers[i], seedLen);
        idxi = (**seeds)[codei].size();

        (**seeds)[codei].push_back(new char[seedLen + 1]);
        strcpy((**seeds)[codei][idxi], kmers[i]);
    }

    if (hitThreshold != -1) {
        for (int i = 0; i < kmersLen; ++i) {

            for (int j = 0; j < seedLen; ++j) {
                char* tmp = new char[seedLen + 1];
                strcpy(tmp, kmers[i]);

                for (int k = 0; k < AA; ++k) {
                    if (AMINO_ACIDS[k] == kmers[i][j]) continue;

                    tmp[j] = AMINO_ACIDS[k];
                    if (permutation(kmers[i], tmp, seedLen, scorer) >= hitThreshold) {
                        codei = seedCode(kmers[i], seedLen);
                        idxi = (**seeds)[codei].size();

                        (**seeds)[codei].push_back(new char[seedLen + 1]);
                        strcpy((**seeds)[codei][idxi], tmp);
                    }
                }

                delete[] tmp;
            }
        }
    }

    kmersDelete(kmers, kmersLen);
}

extern void seedsDelete(Seeds* seeds) {
    for (int i = 0; i < (int) seeds->size(); ++i) {
        for (int j = 0; j < (int) (*seeds)[i].size(); ++j) {
            delete[] (*seeds)[i][j];
        }
        vector<char*>().swap((*seeds)[i]);
    }

    seeds->clear();
    delete seeds;
}

extern void seedScoresCreate(int** seedScores, int* seedScoresLen, int seedLen,
    Scorer* scorer) {

    (*seedScoresLen) = SEED_IDX_LEN(seedLen);
    (*seedScores) = new int[*seedScoresLen];

    vector<char> seed;

    seedScoresCreateRec(seedScores, &seed, seedLen, scorer);
}

extern void seedScoresDelete(int* seedScores) {
    delete[] seedScores;
}

extern void dataCreate(Data** data, int len) {
    vector<int> vi;
    (*data) = new Data(len, vi);
}

extern void dataDelete(Data* data) {
    for (unsigned int i = 0; i < data->size(); ++i) {
        vector<int>().swap((*data)[i]);
    }
    data->clear();
    delete data;
}

// ***************************************************************************

// ***************************************************************************
// PRIVATE

static int seedCode(char* seed, int seedLen) {

    int code = 0;
    int start = 5 * (seedLen - 1);

    for (int i = 0; i < seedLen; ++i) {
        code += static_cast<int>(toupper(seed[i] - 'A')) << (start - 5 * i);
    }

    return code;
}

static int seedCode(vector<char>* seed) {

    int code = 0;
    int start = 5 * (seed->size() - 1);

    for (unsigned int i = 0; i < seed->size(); ++i) {
        code += static_cast<int>((*seed)[i] - 'A') << (start - 5 * i);
    }

    return code;
}

static void kmersCreate(char*** kmers, int* kmersLen, int seedLen) {

    (*kmersLen) = pow(AA, seedLen);
    (*kmers) = new char*[*kmersLen];

    vector<char> kmer;

    kmersCreateRec(kmers, 0, &kmer, seedLen, seedLen);
}

static int kmersCreateRec(char*** kmers, int kmerIdx, vector<char>* kmer,
    int n, int seedLen) {

    if (n == 0) {
        (*kmers)[kmerIdx] = new char[seedLen + 1];
        for (int i = 0; i < seedLen; ++i) {
            (*kmers)[kmerIdx][i] = (*kmer)[i];
        }
        (*kmers)[kmerIdx][seedLen] = '\0';
        return kmerIdx + 1;
    }

    for (int i = 0; i < AA; ++i) {
        kmer->push_back(AMINO_ACIDS[i]);
        kmerIdx = kmersCreateRec(kmers, kmerIdx, kmer, n - 1, seedLen);
        kmer->pop_back();
    }

    return kmerIdx;
}

static void kmersDelete(char** kmers, int kmersLen) {
    for (int i = 0; i < kmersLen; ++i) {
        delete[] kmers[i];
    }
    delete[] kmers;
}

static int permutation(char* kmer1, char* kmer2, int seedLen, Scorer* scorer) {

    int score = 0;

    for (int i = 0; i < seedLen; ++i) {
        score += scorerScore(scorer, scorerEncode(kmer1[i]),
            scorerEncode(kmer2[i]));
    }

    return score;
}

static int seedScore(vector<char>* seed, Scorer* scorer) {

    int score = 0;
    int encodedChar;

    for (int i = 0; i < (int) seed->size(); ++i) {
        encodedChar = scorerEncode((*seed)[i]);
        score += scorerScore(scorer, encodedChar, encodedChar);
    }

    return score;
}

static void seedScoresCreateRec(int** seedScores, vector<char>* seed, int n, Scorer* scorer) {

    if (n == 0) {
        (*seedScores)[seedCode(seed)] = seedScore(seed, scorer);
        return;
    }

    for (int i = 0; i < AA; ++i) {
        (*seed).push_back(AMINO_ACIDS[i]);
        seedScoresCreateRec(seedScores, seed, n - 1, scorer);
        (*seed).pop_back();
    }
}

// ***************************************************************************
