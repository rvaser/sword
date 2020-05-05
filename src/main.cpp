#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <assert.h>
#include <memory>

#include "thread_pool/thread_pool.hpp"

#include "writer.hpp"
#include "evalue.hpp"
#include "score_matrix.hpp"
#include "database_search.hpp"
#include "database_alignment.hpp"
#include "utils.hpp"

static const char* version = "v1.0.4";

static struct option options[] = {
    {"query", required_argument, 0, 'i'},
    {"target", required_argument, 0, 'j'},
    {"gap-extend", required_argument, 0, 'e'},
    {"gap-open", required_argument, 0, 'g'},
    {"matrix", required_argument, 0, 'm'},
    {"out", required_argument, 0, 'o'},
    {"outfmt", required_argument, 0, 'f'},
    {"evalue", required_argument, 0, 'v'},
    {"max-aligns", required_argument, 0, 'a'},
    {"algorithm", required_argument, 0, 'A'},
    {"kmer-length", required_argument, 0, 'k'},
    {"max-candidates", required_argument, 0, 'c'},
    {"threshold", required_argument, 0, 'T'},
    {"threads", required_argument, 0, 't'},
    {"version", no_argument, 0, 'V'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
};

OutputType strToOutputType(const std::string& str);

ScoreMatrixType strToScorerType(const std::string& str);

AlignmentType strToAlignmentType(const std::string& str);

void help();

int main(int argc, char* argv[]) {

    auto threads = std::thread::hardware_concurrency() / 2;

    int32_t gap_open = 10;
    int32_t gap_extend = 1;

    std::string queries_path;
    std::string database_path;

    ScoreMatrixType scorer_type = ScoreMatrixType::kBlosum62;

    std::string output_path;
    OutputType output_format = OutputType::kBm9;

    double max_evalue = 10;
    int32_t max_alignments = 10;

    AlignmentType algorithm = AlignmentType::kSW;

    uint32_t kmer_length = 3;
    uint32_t max_candidates = 30000;
    uint32_t threshold = 13;

    char opt;
    while ((opt = getopt_long(argc, argv, "i:j:g:e:m:o:f:v:a:A:k:c:T:t:h", options, nullptr)) != -1) {

        switch (opt) {
        case 'i':
            queries_path = optarg;
            break;
        case 'j':
            database_path = optarg;
            break;
        case 'g':
            gap_open = atoi(optarg);
            break;
        case 'e':
            gap_extend = atoi(optarg);
            break;
        case 'm':
            scorer_type = strToScorerType(optarg);
            break;
        case 'o':
            output_path = optarg;
            break;
        case 'f':
            output_format = strToOutputType(optarg);
            break;
        case 'v':
            max_evalue = atof(optarg);
            break;
        case 'a':
            max_alignments = atoi(optarg);
            break;
        case 'A':
            algorithm = strToAlignmentType(optarg);
            break;
        case 'k':
            kmer_length = atoi(optarg);
            break;
        case 'c':
            max_candidates = atoi(optarg);
            break;
        case 'T':
            threshold = atoi(optarg);
            break;
        case 't':
            threads = atoi(optarg);
            break;
        case 'V':
            printf("%s\n", version);
            return 0;
        case 'h':
            help();
            return 0;
        default:
            return 1;
        }
    }

    if (queries_path.empty()) {
        fprintf(stderr, "[sword::] error: missing input queries file!\n");
        help();
        return 1;
    }
    if (database_path.empty()) {
        fprintf(stderr, "[sword::] error: missing input database file!\n");
        help();
        return 1;
    }

    std::shared_ptr<thread_pool::ThreadPool> thread_pool = thread_pool::createThreadPool(threads);

    std::shared_ptr<ScoreMatrix> scorer = createScoreMatrix(scorer_type,
        gap_open, gap_extend);

    Timer timer;
    timer.start();

    Indexes indexes;
    auto database_cells = searchDatabase(indexes, database_path, queries_path,
        kmer_length, max_candidates, scorer, threshold, thread_pool);

    timer.stop();
    timer.print("database", "search");
    timer.reset();
    timer.start();

    std::shared_ptr<EValue> evalue_params = createEValue(database_cells, scorer);

    std::vector<AlignmentSet> alignments;
    alignDatabase(alignments, algorithm, database_path, queries_path, indexes,
        max_evalue, evalue_params, max_alignments, scorer, output_path,
        output_format, thread_pool);

    timer.stop();
    timer.print("database", "alignment");

    return 0;
}

OutputType strToOutputType(const std::string& str) {

    if (str.compare("bm0") == 0) {
        return OutputType::kBm0;
    } else if (str.compare("bm8") == 0) {
        return OutputType::kBm8;
    } else if (str.compare("bm9") == 0) {
        return OutputType::kBm9;
    }

    assert(false && "unrecognized output format");
}

ScoreMatrixType strToScorerType(const std::string& str) {

    if (str.compare("BLOSUM_45") == 0) {
        return ScoreMatrixType::kBlosum45;
    } else if (str.compare("BLOSUM_50") == 0) {
        return ScoreMatrixType::kBlosum50;
    } else if (str.compare("BLOSUM_62") == 0) {
        return ScoreMatrixType::kBlosum62;
    } else if (str.compare("BLOSUM_80") == 0) {
        return ScoreMatrixType::kBlosum80;
    } else if (str.compare("BLOSUM_90") == 0) {
        return ScoreMatrixType::kBlosum90;
    } else if (str.compare("PAM_30") == 0) {
        return ScoreMatrixType::kPam30;
    } else if (str.compare("PAM_70") == 0) {
        return ScoreMatrixType::kPam70;
    } else if (str.compare("PAM_250") == 0) {
        return ScoreMatrixType::kPam250;
    }

    assert(false && "unrecognized scorer type");
}

AlignmentType strToAlignmentType(const std::string& str) {

    if (str.compare("NW") == 0) {
        return AlignmentType::kNW;
    } else if (str.compare("HW") == 0) {
        return AlignmentType::kHW;
    } else if (str.compare("OV") == 0) {
        return AlignmentType::kOV;
    } else if (str.compare("SW") == 0) {
        return AlignmentType::kSW;
    }

    assert(false && "unrecognized aignment type");
}

void help() {
    printf(
    "usage: sword -i <query db file> -j <target db file> [arguments ...]\n"
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
    "    -m, --matrix <string>\n"
    "        default: BLOSUM_62\n"
    "        similarity matrix, can be one of the following:\n"
    "            BLOSUM_45\n"
    "            BLOSUM_50\n"
    "            BLOSUM_62\n"
    "            BLOSUM_80\n"
    "            BLOSUM_90\n"
    "            PAM_30\n"
    "            PAM_70\n"
    "            PAM_250\n"
    "    -o, --out <string>\n"
    "        default: stdout\n"
    "        output file for the alignment\n"
    "    -f, --outfmt <string>\n"
    "        default: bm9\n"
    "        out format for the output file, must be one of the following:\n"
    "            bm0      - blast m0 output format\n"
    "            bm8      - blast m8 tabular output format\n"
    "            bm9      - blast m9 commented tabular output format\n"
    "    -v, --evalue <float>\n"
    "        default: 10.0\n"
    "        evalue threshold, alignments with higher evalue are filtered,\n"
    "        must be given as a positive float\n"
    "    -a, --max-aligns <int>\n"
    "        default: 10\n"
    "        maximum number of alignments to be outputted\n"
    "    -A, --algorithm <string>\n"
    "        default: SW\n"
    "        algorithm used for alignment, must be one of the following: \n"
    "            SW - Smith-Waterman local alignment\n"
    "            NW - Needleman-Wunsch global alignment\n"
    "            HW - semiglobal alignment\n"
    "            OV - overlap alignment\n"
    "    -k, --kmer-length <int>\n"
    "        default: 3\n"
    "        length of kmers used for database search\n"
    "        possible values: 3, 4, 5\n"
    "    -c, --max-candidates <int>\n"
    "        default: 30000\n"
    "        number of target sequences (per query sequence) passed\n"
    "        to the alignment part\n"
    "    -T, --threshold <int>\n"
    "        default: 13\n"
    "        minimum score for two kmers to trigger a hit\n"
    "        if 0 given, only exact matching kmers are checked\n"
    "    -t, --threads <int>\n"
    "        default: hardware concurrency / 2\n"
    "        number of threads used in thread pool\n"
    "    --version\n"
    "        prints the version number\n"
    "    -h, --help\n"
    "        prints out the help\n");
}
