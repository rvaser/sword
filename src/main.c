#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "timer.h"
#include "database_heuristics.h"
#include "swsharp/evalue.h"
#include "swsharp/swsharp.h"

#define ASSERT(expr, fmt, ...)\
    do {\
        if (!(expr)) {\
            fprintf(stderr, "[ERROR]: " fmt "\n", ##__VA_ARGS__);\
            exit(-1);\
        }\
    } while (0)

#define CHAR_INT_LEN(x) (sizeof(x) / sizeof(CharInt))

typedef struct CharInt {
    const char* format;
    const int code;
} CharInt;

typedef struct ValueFunctionParam {
    Scorer* scorer;
    int totalLength;
} ValueFunctionParam;

static struct option options[] = {
    {"cards", required_argument, 0, 'c'},
    {"gap-extend", required_argument, 0, 'e'},
    {"gap-open", required_argument, 0, 'g'},
    {"query", required_argument, 0, 'i'},
    {"target", required_argument, 0, 'j'},
    {"matrix", required_argument, 0, 'm'},
    {"out", required_argument, 0, 'o'},
    {"outfmt", required_argument, 0, 'f'},
    {"evalue", required_argument, 0, 'E'},
    {"max-aligns", required_argument, 0, 'M'},
    {"algorithm", required_argument, 0, 'A'},
    {"nocache", no_argument, 0, 'C'},
    {"cpu", no_argument, 0, 'P'},
    {"threads", required_argument, 0, 'T'},
    {"seed-length", required_argument, 0, 's'},
    {"hit-threshold", required_argument, 0, 't'},
    {"max-candidates", required_argument, 0, 'd'},
    {"dbchunk", required_argument, 0, 'D'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
};

static CharInt outFormats[] = {
    { "bm0", SW_OUT_DB_BLASTM0 },
    { "bm8", SW_OUT_DB_BLASTM8 },
    { "bm9", SW_OUT_DB_BLASTM9 },
    { "light", SW_OUT_DB_LIGHT }
};

static CharInt algorithms[] = {
    { "SW", SW_ALIGN },
    { "NW", NW_ALIGN },
    { "HW", HW_ALIGN },
    { "OV", OV_ALIGN }
};

static void help();

static void getCudaCards(int** cards, int* cardsLen, char* optarg);

static int getOutFormat(char* optarg);
static int getAlgorithm(char* optarg);

static void valueFunction(double* values, int* scores, Chain* query,
    Chain** database, int databaseLen, int* cards, int cardsLen, void* param);

int main(int argc, char* argv[]) {

    char* queryPath = NULL;
    char* databasePath = NULL;

    int gapOpen = 10;
    int gapExtend = 1;

    char* matrix = "BLOSUM_62";

    int maxAlignments = 10;
    float maxEValue = 10;

    int cardsLen = -1;
    int* cards = NULL;

    char* out = NULL;
    int outFormat = SW_OUT_DB_BLASTM9;

    int algorithm = SW_ALIGN;

    int cache = 1;

    int forceCpu = 0;

    int threads = 8;

    int seedLen = 3;
    int hitThreshold = 13;

    int maxCandidates = 30000;

    float dbchunk = 2;

    while (1) {

        char argument = getopt_long(argc, argv, "i:j:g:e:s:t:T:h", options, NULL);

        if (argument == -1) {
            break;
        }

        switch (argument) {
        case 'i':
            queryPath = optarg;
            break;
        case 'j':
            databasePath = optarg;
            break;
        case 'g':
            gapOpen = atoi(optarg);
            break;
        case 'e':
            gapExtend = atoi(optarg);
            break;
        case 'c':
            getCudaCards(&cards, &cardsLen, optarg);
            break;
        case 'o':
            out = optarg;
            break;
        case 'f':
            outFormat = getOutFormat(optarg);
            break;
        case 'M':
            maxAlignments = atoi(optarg);
            break;
        case 'E':
            maxEValue = atof(optarg);
            break;
        case 'm':
            matrix = optarg;
            break;
        case 'A':
            algorithm = getAlgorithm(optarg);
            break;
        case 'C':
            cache = 0;
            break;
        case 'P':
            forceCpu = 1;
            break;
        case 'T':
            threads = atoi(optarg);
            break;
        case 's':
            seedLen = atoi(optarg);
            break;
        case 't':
            hitThreshold = atoi(optarg);
            break;
        case 'd':
            maxCandidates = atoi(optarg);
            break;
        case 'D':
            dbchunk = atof(optarg);
            break;
        default:
            help();
            return -1;
        }
    }

    ASSERT(queryPath != NULL, "missing option -i (query file)");
    ASSERT(databasePath != NULL, "missing option -j (database file)");

    ASSERT(seedLen > 2 && seedLen < 6, "seed length possible values = 3,4,5");
    ASSERT(maxCandidates > 0, "max-candidates possible values = [1, infinity>");

    ASSERT(dbchunk > 0, "invalid dbchunk size");

    if (forceCpu) {
        cards = NULL;
        cardsLen = 0;
    } else {

        if (cardsLen == -1) {
            cudaGetCards(&cards, &cardsLen);
        }

        ASSERT(cudaCheckCards(cards, cardsLen), "invalid cuda cards");
    }

    ASSERT(maxEValue > 0, "invalid evalue");

    ASSERT(threads >= 0, "invalid thread number");
    threadPoolInitialize(threads);

    Scorer* scorer = NULL;
    scorerCreateMatrix(&scorer, matrix, gapOpen, gapExtend);

    void* indices = databaseIndicesCreate(databasePath, queryPath, seedLen, maxCandidates,
        hitThreshold, scorer, threads);

    Timeval swTimer;
    timerStart(&swTimer);

    if (indices != NULL) {

        Chain** queries = NULL;
        int queriesLen = 0;
        readFastaChains(&queries, &queriesLen, queryPath);

        if (cache) {
            dumpFastaChains(databasePath);
        }

        int chains;
        long long cells;
        statFastaChains(&chains, &cells, databasePath);

        EValueParams* eValueParams = createEValueParams(cells, scorer);

        DbAlignment*** dbAlignments = NULL;
        int* dbAlignmentsLens = NULL;

        int dbReadSize = (int) (dbchunk * 1000000000); // * ~1GB

        Chain** database = NULL;
        int databaseLen = 0;
        int databaseStart = 0;

        FILE* handle;
        int serialized;

        readFastaChainsPartInit(&database, &databaseLen, &handle, &serialized, databasePath);

        int i, j;

        while (1) {

            int status = readFastaChainsPart(&database, &databaseLen, handle,
                serialized, dbReadSize);

            DbAlignment*** dbAlignmentsPart =
                (DbAlignment***) malloc(queriesLen * sizeof(DbAlignment**));
            int* dbAlignmentsPartLens = (int*) malloc(queriesLen * sizeof(int));

            Chain** filteredDatabase = NULL;
            int filteredDatabaseLen = 0;

            ChainDatabase* chainDatabase = NULL;

            int* usedIndices = NULL;
            char* usedMask = (char*) calloc(databaseLen, sizeof(char));

            for (i = 0; i < queriesLen; ++i) {

                usedIndices = filteredDatabaseCreate(&filteredDatabase,
                    &filteredDatabaseLen, indices, i, database, databaseLen,  1);

                if (filteredDatabaseLen == 0) {
                    dbAlignmentsPart[i] = NULL;
                    dbAlignmentsPartLens[i] = 0;
                    continue;
                }

                chainDatabase = chainDatabaseCreate(filteredDatabase, 0,
                    filteredDatabaseLen, cards, cardsLen);

                alignDatabase(&dbAlignmentsPart[i], &dbAlignmentsPartLens[i], algorithm, queries[i],
                    chainDatabase, scorer, maxAlignments, valueFunction, (void*) eValueParams,
                    maxEValue, NULL, 0, cards, cardsLen, NULL);

                if (usedIndices != NULL) {
                    for (j = 0; j < dbAlignmentsPartLens[i]; ++j) {

                        DbAlignment* dbAlignment = dbAlignmentsPart[i][j];
                        int targetIdx = dbAlignmentGetTargetIdx(dbAlignment);

                        usedMask[usedIndices[targetIdx]] = 1;
                    }

                    free(usedIndices);
                }

                chainDatabaseDelete(chainDatabase);

                filteredDatabaseDelete(filteredDatabase);
            }

            if (dbAlignments == NULL) {
                dbAlignments = dbAlignmentsPart;
                dbAlignmentsLens = dbAlignmentsPartLens;
            } else {
                dbAlignmentsMerge(dbAlignments, dbAlignmentsLens, dbAlignmentsPart,
                    dbAlignmentsPartLens, queriesLen, maxAlignments);
                deleteShotgunDatabase(dbAlignmentsPart, dbAlignmentsPartLens, queriesLen);
            }

            for (i = databaseStart; i < databaseLen; ++i) {
                if (!usedMask[i] && database[i] != NULL) {
                    chainDelete(database[i]);
                    database[i] = NULL;
                }
            }

            free(usedMask);

            if (status == 0) {
                break;
            }

            databaseStart = databaseLen;
        }

        fclose(handle);

        outputShotgunDatabase(dbAlignments, dbAlignmentsLens, queriesLen, out, outFormat);
        deleteShotgunDatabase(dbAlignments, dbAlignmentsLens, queriesLen);

        deleteEValueParams(eValueParams);

        deleteFastaChains(database, databaseLen);

        deleteFastaChains(queries, queriesLen);

        databaseIndicesDelete(indices);
    }

    timerPrint("Sw#Part", timerStop(&swTimer));

    scorerDelete(scorer);

    threadPoolTerminate();
    free(cards);

    return 0;
}

static void getCudaCards(int** cards, int* cardsLen, char* optarg) {

    *cardsLen = strlen(optarg);
    *cards = (int*) malloc(*cardsLen * sizeof(int));

    int i;
    for (i = 0; i < *cardsLen; ++i) {
        (*cards)[i] = optarg[i] - '0';
    }
}

static int getOutFormat(char* optarg) {

    int i;
    for (i = 0; i < CHAR_INT_LEN(outFormats); ++i) {
        if (strcmp(outFormats[i].format, optarg) == 0) {
            return outFormats[i].code;
        }
    }

    ASSERT(0, "unknown out format %s", optarg);
}

static int getAlgorithm(char* optarg) {

    int i;
    for (i = 0; i < CHAR_INT_LEN(algorithms); ++i) {
        if (strcmp(algorithms[i].format, optarg) == 0) {
            return algorithms[i].code;
        }
    }

    ASSERT(0, "unknown algorithm %s", optarg);
}

static void valueFunction(double* values, int* scores, Chain* query,
    Chain** database, int databaseLen, int* cards, int cardsLen, void* param_ ) {

    EValueParams* eValueParams = (EValueParams*) param_;
    eValues(values, scores, query, database, databaseLen, cards, cardsLen, eValueParams);
}

static void help() {
    printf(
    "usage: swsharpdbh -i <query db file> -j <target db file> [arguments ...]\n"
    "\n"
    "arguments:\n"
    "    -i, --query <file>\n"
    "        (required)\n"
    "        input fasta database query file\n"
    "    -j, --target <file>\n"
    "        (required)\n"
    "        input fasta database target file\n"
    "    -g, --gap-open <int>\n"
    "        default: 10\n"
    "        gap opening penalty, must be given as a positive integer \n"
    "    -e, --gap-extend <int>\n"
    "        default: 1\n"
    "        gap extension penalty, must be given as a positive integer and\n"
    "        must be less or equal to gap opening penalty\n"
    "    --matrix <string>\n"
    "        default: BLOSUM_62\n"
    "        similarity matrix, can be one of the following:\n"
    "            BLOSUM_45\n"
    "            BLOSUM_50\n"
    "            BLOSUM_62\n"
    "            BLOSUM_80\n"
    "            BLOSUM_90\n"
    "            BLOSUM_30\n"
    "            BLOSUM_70\n"
    "            BLOSUM_250\n"
    "            EDNA_FULL\n"
    "    --evalue <float>\n"
    "        default: 10.0\n"
    "        evalue threshold, alignments with higher evalue are filtered,\n"
    "        must be given as a positive float\n"
    "    --max-aligns <int>\n"
    "        default: 10\n"
    "        maximum number of alignments to be outputted\n"
    "    --algorithm <string>\n"
    "        default: SW\n"
    "        algorithm used for alignment, must be one of the following: \n"
    "            SW - Smith-Waterman local alignment\n"
    "            NW - Needleman-Wunsch global alignment\n"
    "            HW - semiglobal alignment\n"
    "            OV - overlap alignment\n"
    "    --cards <ints>\n"
    "        default: all available CUDA cards\n"
    "        list of cards should be given as an array of card indexes delimited with\n"
    "        nothing, for example usage of first two cards is given as --cards 01\n"
    "    --out <string>\n"
    "        default: stdout\n"
    "        output file for the alignment\n"
    "    --outfmt <string>\n"
    "        default: bm9\n"
    "        out format for the output file, must be one of the following:\n"
    "            bm0      - blast m0 output format\n"
    "            bm8      - blast m8 tabular output format\n"
    "            bm9      - blast m9 commented tabular output format\n"
    "            light    - score-name tabbed output\n"
    "    --nocache\n"
    "        serialized database is stored to speed up future runs with the\n"
    "        same database, option disables this behaviour\n"
    "    --cpu\n"
    "        only cpu is used\n"
    "    -T --threads <int>\n"
    "        default: 8\n"
    "        number of threads used in thread pool\n"
    "    -s, --seed-length <int>\n"
    "        default: 3\n"
    "        length of seeds for hash search\n"
    "        possible values are: 3,4,5\n"
    "    -t, --hit-threshold <int>\n"
    "        default: 13\n"
    "        minimum score for two seeds to trigger a hit\n"
    "        use -1 to check only exact matching seeds\n"
    "    --max-candidates <int>\n"
    "        default: 30000\n"
    "        number of database sequences passed on to the Smith-Waterman part\n"
    "    --dbchunk <float>\n"
    "        default: 2\n"
    "        size of database (in GB) to be read and aligned at a time"
    "    -h, -help\n"
    "        prints out the help\n");
}
