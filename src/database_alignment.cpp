/*!
 * @file database_alignment.cpp
 *
 * @brief Database alignment source file
 */

#include <assert.h>
#include <algorithm>

#include "opal.h"

#include "chain.hpp"
#include "reader.hpp"
#include "writer.hpp"
#include "score_matrix.hpp"
#include "evalue.hpp"
#include "database_alignment.hpp"

constexpr size_t kDatabasePartSize = 1000000000; /* ~1 GB */

std::unique_ptr<Alignment> createAlignment(int32_t score, double evalue,
    uint32_t query_id, uint32_t target_id) {

    return std::unique_ptr<Alignment>(new Alignment(score, evalue, query_id,
        target_id));
}

Alignment::Alignment(int32_t score, double evalue, uint32_t query_id, uint32_t target_id)
        : score_(score), evalue_(evalue), query_id_(query_id), target_id_(target_id) {
}

void Alignment::update(uint32_t query_begin, uint32_t query_end, uint32_t target_begin,
    uint32_t target_end, const unsigned char* alignment, uint32_t length) {

    query_begin_ = query_begin;
    query_end_ = query_end;
    target_begin_ = target_begin;
    target_end_ = target_end;

    alignment_.clear();
    for (uint32_t i = 0; i < length; ++i) {
        alignment_ += alignment[i];
    }
}

/* ************************************************************************** */
/* Preprocess */

bool compareAlignment(const std::unique_ptr<Alignment>& left,
    const std::unique_ptr<Alignment>& right) {

    if (left->evalue() < right->evalue()) return true;
    if (left->evalue() == right->evalue() && left->score() > right->score()) return true;
    return false;
}

unsigned char* strToUnsignedCharPtr(const std::string& src) {
    auto dst = new unsigned char[src.size()];
    for (uint32_t i = 0; i < src.size(); ++i) {
        dst[i] = src[i];
    }
    return dst;
}

uint32_t alignmentTypeToOpalMode(AlignmentType algorithm) {

    switch (algorithm) {
        case AlignmentType::kNW:
            return OPAL_MODE_NW;
        case AlignmentType::kHW:
            return OPAL_MODE_HW;
        case AlignmentType::kOV:
            return OPAL_MODE_OV;
        case AlignmentType::kSW:
        default:
            return OPAL_MODE_SW;
    }
}

/* ************************************************************************** */

void scoreChains(AlignmentSet& dst, const std::unique_ptr<Chain>& query,
    std::vector<uint32_t>& indexes, const ChainSet& database, uint32_t database_start,
    uint32_t algorithm, double max_evalue, std::shared_ptr<EValue> evalue_params,
    uint32_t max_alignments, std::shared_ptr<ScoreMatrix> scorer) {

    uint32_t i;
    for (i = 0; i < indexes.size(); ++i) {
        if (indexes[i] >= database.size()) break;
    }

    auto database_length = i;
    if (database_length == 0) return;

    OpalSearchResult* results[database_length];
    for (i = 0; i < database_length; ++i) {
        results[i] = new OpalSearchResult();
        opalInitSearchResult(results[i]);
    }

    unsigned char* query_ = strToUnsignedCharPtr(query->data());
    int query_length = query->length();

    unsigned char* database_[database_length];
    int database_lengths[database_length];

    for (i = 0; i < database_length; ++i) {
        const auto& target = database[indexes[i]];
        database_[i] = strToUnsignedCharPtr(target->data());
        database_lengths[i] = target->length();
    }

    auto error = opalSearchDatabase(query_, query_length, database_, database_length,
        database_lengths, scorer->gap_open(), scorer->gap_extend(), scorer->data(),
        scorer->num_rows_, results, OPAL_SEARCH_SCORE, algorithm,
        OPAL_OVERFLOW_SIMPLE);

    if (error) {
        fprintf(stderr, "Opal alignment failed with code %d\n", error);
    }

    for (i = 0; i < database_length; ++i) {
        if (results[i]->scoreSet == 1) {

            auto evalue = evalue_params->calculate(results[i]->score,
                query_length, database_lengths[i]);

            if (evalue <= max_evalue) {
                dst.emplace_back(createAlignment(results[i]->score, evalue,
                    query->id(), indexes[i]));
            }
        }
    }

    std::sort(dst.begin(), dst.end(), compareAlignment);

    if (max_alignments && dst.size() > max_alignments) {
        dst.resize(max_alignments);
    }

    std::vector<uint32_t> temp(indexes.begin() + database_length, indexes.end());
    indexes.swap(temp);

    for (const auto& it: database_) {
        delete[] it;
    }

    delete[] query_;

    for (auto& it: results) {
        delete it;
    }
}

void alignChains(AlignmentSet& dst, const std::unique_ptr<Chain>& query,
    const ChainSet& database, uint32_t algorithm,
    std::shared_ptr<ScoreMatrix> scorer) {

    if (dst.size() == 0) return;

    auto database_length = dst.size();

    OpalSearchResult* results[database_length];
    for (uint32_t i = 0; i < database_length; ++i) {
        results[i] = new OpalSearchResult();
        opalInitSearchResult(results[i]);
    }

    unsigned char* query_ = strToUnsignedCharPtr(query->data());
    int query_length = query->length();

    unsigned char* database_[database_length];
    int database_lengths[database_length];

    for (uint32_t i = 0; i < database_length; ++i) {
        const auto& target = database[dst[i]->target_id()];
        database_[i] = strToUnsignedCharPtr(target->data());
        database_lengths[i] = target->length();
    }

    auto error = opalSearchDatabase(query_, query_length, database_, database_length,
        database_lengths, scorer->gap_open(), scorer->gap_extend(), scorer->data(),
        scorer->num_rows_, results, OPAL_SEARCH_ALIGNMENT,
        algorithm, OPAL_OVERFLOW_SIMPLE);

    if (error) {
        fprintf(stderr, "Opal alignment failed with code %d\n", error);
    }

    for (uint32_t i = 0; i < database_length; ++i) {
        dst[i]->update(results[i]->startLocationQuery, results[i]->endLocationQuery,
            results[i]->startLocationTarget, results[i]->endLocationTarget,
            results[i]->alignment, results[i]->alignmentLength);
    }

    for (const auto& it: database_) {
        delete[] it;
    }

    delete[] query_;

    for (auto& it: results) {
        free(it->alignment);
        delete it;
    }
}

void alignDatabase(std::vector<AlignmentSet>& dst, AlignmentType algorithm_,
    const std::string& database_path, const std::string& queries_path,
    Indexes& indexes, double max_evalue, std::shared_ptr<EValue> evalue_params,
    uint32_t max_alignments, std::shared_ptr<ScoreMatrix> scorer,
    const std::string& output_path, OutputType output_format,
    std::shared_ptr<thread_pool::ThreadPool> thread_pool) {

    auto algorithm = alignmentTypeToOpalMode(algorithm_);

    ChainSet queries;
    createChainSet(queries, queries_path);

    dst.clear();
    dst.resize(queries.size());

    ChainSet database;
    uint32_t database_start = 0;
    std::shared_ptr<Reader> reader = createChainSetPartInitialize(database_path);

    /* find scores for indexed targets */
    while (true) {

        auto status = createChainSetPart(database, reader, kDatabasePartSize);

        std::vector<std::future<void>> thread_futures;

        for (uint32_t i = 0; i < queries.size(); ++i) {
            thread_futures.emplace_back(thread_pool->submit(scoreChains,
                std::ref(dst[i]), std::ref(queries[i]), std::ref(indexes[i]),
                std::ref(database), database_start, algorithm, max_evalue,
                evalue_params, max_alignments, scorer));
        }

        for (const auto& it: thread_futures) {
            it.wait();
        }

        auto used_mask = new uint8_t[database.size()]();
        for (const auto& it: dst) {
            for (const auto& alignment: it) {
                used_mask[alignment->target_id()] = 1;
            }
        }

        for (uint32_t i = database_start; i < database.size(); ++i) {
            if (used_mask[i] == 0) {
                database[i].reset(nullptr);
            }
        }

        delete[] used_mask;

        if (status == false) {
            break;
        }

        database_start = database.size();
    }

    /* find alignments for best targets */
    {
        std::vector<std::future<void>> thread_futures;

        for (uint32_t i = 0; i < queries.size(); ++i) {
            thread_futures.emplace_back(thread_pool->submit(alignChains,
                std::ref(dst[i]), std::ref(queries[i]), std::ref(database),
                algorithm, scorer));
        }

        for (const auto& it: thread_futures) {
            it.wait();
        }
    }

    auto writer = createWriter(output_path, output_format, scorer);

    for (const auto& it: dst) {
        writer->write_alignments(it, queries, database);
    }
}
