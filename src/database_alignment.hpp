/*!
 * @file database_alignment.hpp
 *
 * @brief Database alignment header file
 */

#pragma once

#include <stdint.h>
#include <string.h>
#include <vector>
#include <memory>

#include "thread_pool/thread_pool.hpp"

enum class OutputType;
class ScoreMatrix;
class EValue;
class Alignment;

enum class AlignmentType {
    kNW, // global alignment (Needleman-Wunsch)
    kHW, // semi-global alignemnt (query - _target_)
    kOV, // semi-global alignment (_query_ - _target_)
    kSW // local alignment (Smith-Waterman)
};

using Indexes = std::vector<std::vector<uint32_t>>;
using AlignmentSet = std::vector<std::unique_ptr<Alignment>>;

std::unique_ptr<Alignment> createAlignment(int32_t score, double evalue,
    uint32_t query_id, uint32_t target_id);

class Alignment {
public:

    ~Alignment() = default;

    int32_t score() const {
        return score_;
    }

    double evalue() const {
        return evalue_;
    }

    uint32_t query_id() const {
        return query_id_;
    }

    uint32_t query_begin() const {
        return query_begin_;
    }

    uint32_t query_end() const {
        return query_end_;
    }

    uint32_t target_id() const {
        return target_id_;
    }

    uint32_t target_begin() const {
        return target_begin_;
    }

    uint32_t target_end() const {
        return target_end_;
    }

    const std::string& alignment() const {
        return alignment_;
    }

    void update(uint32_t query_begin, uint32_t query_end, uint32_t target_begin,
        uint32_t target_end, const unsigned char* alignment, uint32_t length);

    friend std::unique_ptr<Alignment> createAlignment(int32_t score,
        double evalue, uint32_t query_id, uint32_t target_id);

private:

    Alignment(int32_t score, double evalue, uint32_t query_id, uint32_t target_id);

	Alignment(const Alignment&) = delete;
	const Alignment& operator=(const Alignment&) = delete;

    int32_t score_;
    double evalue_;

    uint32_t query_id_;
    uint32_t target_id_;

    uint32_t query_begin_;
    uint32_t query_end_;
    uint32_t target_begin_;
    uint32_t target_end_;

    std::string alignment_;
};

void alignDatabase(std::vector<AlignmentSet>& dst, AlignmentType algorithm,
    const std::string& database_path, const std::string& queries_path,
    Indexes& indexes, double max_evalue, std::shared_ptr<EValue> evalue_params,
    uint32_t max_alignments, std::shared_ptr<ScoreMatrix> scorer,
    const std::string& output_path, OutputType output_format,
    std::shared_ptr<thread_pool::ThreadPool> thread_pool);
