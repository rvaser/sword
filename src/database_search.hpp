/*!
 * @file database_search.hpp
 *
 * @brief Database search header file
 */

#pragma once

#include <stdint.h>
#include <vector>
#include <string>
#include <memory>

#include "thread_pool/thread_pool.hpp"

class ScoreMatrix;

using Indexes = std::vector<std::vector<uint32_t>>;

uint64_t searchDatabase(Indexes& dst, const std::string& database_path,
    const std::string& queries_path, uint32_t kmer_length, uint32_t max_candidates,
    std::shared_ptr<ScoreMatrix> score_matrix, uint32_t score_threshold,
    std::shared_ptr<thread_pool::ThreadPool> thread_pool);
