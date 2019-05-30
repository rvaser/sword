/*!
 * @file database_search.cpp
 *
 * @brief Database search source file
 */

#include <algorithm>

#include "chain.hpp"
#include "reader.hpp"
#include "kmers.hpp"
#include "hash.hpp"
#include "utils.hpp"
#include "database_search.hpp"

constexpr size_t kDatabasePartSize = 1000000000; /* ~1 GB */
constexpr uint32_t kMaxShortChainLength = 2000;

constexpr uint32_t kProtBits = 5;
std::vector<uint32_t> kKmerDelMask = { 0, 0, 0, 0x7FFF, 0xFFFFF, 0x1FFFFFF };

/* ************************************************************************** */
/* ChainEntry - used to store additional data of Chain objects */

class ChainEntry;
using ChainEntrySet = std::vector<std::vector<ChainEntry>>;

class ChainEntry {
public:

    ChainEntry() {};
    ChainEntry(uint32_t chain_idx, uint32_t data)
            : chain_idx_(chain_idx), data_(data) {
    };

    ~ChainEntry() = default;

    uint32_t chain_idx() const {
        return chain_idx_;
    }

    uint32_t data() const {
        return data_;
    }

private:

    uint32_t chain_idx_;
    uint32_t data_;
};

bool compareChainEntryDsc(const ChainEntry& left, const ChainEntry& right) {
    return left.data() > right.data();
}

/* ************************************************************************** */
/* Chain preproces */

bool compareChainByLengthAsc(const std::unique_ptr<Chain>& left,
    const std::unique_ptr<Chain>& right) {

    return left->length() < right->length();
}

void preprocDatabase(std::vector<uint32_t>& dst, ChainSet& database,
    size_t num_threads) {

    sort(database.begin(), database.end(), compareChainByLengthAsc);

    uint64_t short_total_length = 0;
    uint64_t long_total_length = 0;
    uint32_t split = 0;

    for (uint32_t i = 0; i < database.size(); ++i) {
        if (database[i]->length() > kMaxShortChainLength) {
            if (split == 0) split = i;
            long_total_length += database[i]->length();
        } else {
            short_total_length += database[i]->length();
        }
    }

    if (short_total_length == 0) {
        split = 0;
    }

    if (long_total_length == 0) {
        split = database.size();
    }

    uint64_t short_task_size = short_total_length / num_threads;
    uint64_t long_task_size = long_total_length / num_threads;

    dst.reserve(2 * num_threads + 1);
    dst.emplace_back(0);

    uint64_t total_length = 0;

    for (uint32_t i = 0; i < split; ++i) {
        total_length += database[i]->length();

        if (total_length > short_task_size) {
            total_length = 0;
            dst.emplace_back(i + 1);
        }
    }

    if (dst.back() != split) {
        dst.emplace_back(split);
    }

    total_length = 0;
    for (uint32_t i = split; i < database.size(); ++i) {
        total_length += database[i]->length();

        if (total_length > long_task_size) {
            total_length = 0;
            dst.emplace_back(i + 1);
        }
    }

    if (dst.back() != database.size()) {
        dst.emplace_back(database.size());
    }
}

/* ************************************************************************** */

using MutexPtr = std::unique_ptr<std::mutex>;

void scoreChains(ChainEntrySet& dst, std::vector<MutexPtr>& entry_mutexes,
    size_t max_candidates, const ChainSet& queries, const ChainSet& database,
    uint32_t database_start, uint32_t database_end, std::shared_ptr<Kmers> kmers) {

    timeval start;
    gettimeofday(&start, nullptr);

    ChainEntrySet entries_part(queries.size());

    std::unique_ptr<uint16_t[]> min_entry_score(new uint16_t[queries.size()]);
    std::unique_ptr<uint16_t[]> entries_found(new uint16_t[queries.size()]);

    {
        for (uint32_t i = 0; i < queries.size(); ++i) {
            auto id = queries[i]->id();
            std::unique_lock<std::mutex> lock(*(entry_mutexes[id].get()));
            entries_found[i] = dst[id].size();
            min_entry_score[i] = entries_found[i] > 0 ? dst[id].back().data() : 65000;
        }
    }

    uint32_t kmer_length = kmers->kmer_length();
    uint32_t max_target_length = database[database_end - 1]->length();
    size_t groups = 0;

    uint32_t kmer_offset = kmer_length - 1;
    uint32_t del_mask = kKmerDelMask[kmer_length];

    uint32_t max_scores_length = kmer_length == 3 ? 100000 : 500000;
    std::unique_ptr<uint16_t[]> scores(new uint16_t[max_scores_length]());
    std::unique_ptr<uint32_t[]> score_lengths(new uint32_t[queries.size()]);
    std::unique_ptr<uint32_t[]> score_starts(new uint32_t[queries.size()+1]);
    score_starts[0] = 0;
    std::unique_ptr<uint16_t[]> max_score(new uint16_t[queries.size()]());

    uint32_t min_score = kmer_length == 3 ? 1 : 0;

    for (uint32_t i = 0; i < queries.size();) {

        ++groups;

        uint32_t group_length = 0;
        uint32_t scores_length = 0;

        for (uint32_t j = i; j < queries.size(); ++j) {

            uint32_t length = queries[j]->length() + max_target_length -
                2 * kmer_length + 1;

            if (scores_length + length > max_scores_length && group_length > 0) {
                break;
            }

            scores_length += length;
            ++group_length;
        }

        auto hash = createHash(queries, i, group_length, kmers);
        Hash::Iterator begin, end;

        for (uint32_t j = database_start; j < database_end; ++j) {

            for (uint32_t k = 0; k < group_length; ++k) {
                score_lengths[k] = queries[i + k]->length() + database[j]->length()
                    - 2 * kmer_length + 1;
                score_starts[k + 1] = score_starts[k] + score_lengths[k];
            }

            const auto& sequence = database[j]->data();
            uint32_t kmer = sequence[0];
            for (uint32_t k = 1; k < kmer_offset; ++k) {
                kmer = (kmer << kProtBits) | sequence[k];
            }

            uint32_t max_diag_id = database[j]->length() - kmer_length;
            for (uint32_t k = kmer_offset; k < sequence.size(); ++k) {
                kmer = ((kmer << kProtBits) | sequence[k]) & del_mask;
                hash->hits(begin, end, kmer);
                for (; begin != end; ++begin) {
                    auto diagonal = max_diag_id + kmer_offset - k +
                        begin->position() + score_starts[begin->id()];
                    ++scores[diagonal];
                    if (max_score[begin->id()] < scores[diagonal]) {
                        max_score[begin->id()] = scores[diagonal];
                    }
                }
            }

            for (uint32_t k = 0; k < group_length; ++k) {

                if (max_score[k] <= min_score) {
                    continue;
                }

                auto id = queries[i + k]->id();
                auto flag = entries_part[id].size() < max_candidates &&
                    entries_found[k] < max_candidates;

                if (flag || max_score[k] >= min_entry_score[k]) {

                    entries_part[id].emplace_back(database[j]->id(), max_score[k]);

                    if (min_entry_score[k] > max_score[k]) {
                        min_entry_score[k] = max_score[k];
                    }
                }
            }

            for (uint32_t k = 0; k < group_length; ++k) {

                if (max_score[k] == 0) {
                    continue;
                }

                max_score[k] = 0;
                std::fill_n(&scores[0] + score_starts[k], score_lengths[k], 0);
            }
        }

        {
            for (uint32_t k = 0; k < group_length; ++k) {
                auto id = queries[i + k]->id();
                std::unique_lock<std::mutex> lock(*(entry_mutexes[id].get()));

                dst[id].insert(dst[id].end(), entries_part[id].begin(),
                    entries_part[id].end());
                std::vector<ChainEntry>().swap(entries_part[id]);

                std::stable_sort(dst[id].begin(), dst[id].end(), compareChainEntryDsc);

                if (dst[id].size() > max_candidates) {
                    std::vector<ChainEntry> temp(dst[id].begin(), dst[id].begin() + max_candidates);
                    dst[id].swap(temp);
                }
            }
        }

        i += group_length;
    }

    timeval stop;
    gettimeofday(&stop, nullptr);

    uint64_t time_ = ((stop.tv_sec - start.tv_sec) * 1000000L + stop.tv_usec) - start.tv_usec;
    fprintf(stderr, "[%u]-[%u]: max = %u, groups = %zu, time = %.5lf s\n",
        database_start, database_end, max_target_length, groups, time_ / (double) 1000000);
}

uint64_t searchDatabase(Indexes& dst, const std::string& database_path,
    const std::string& queries_path, uint32_t kmer_length, uint32_t max_candidates,
    std::shared_ptr<ScoreMatrix> score_matrix, uint32_t score_threshold,
    std::shared_ptr<thread_pool::ThreadPool> thread_pool) {

    ChainSet queries;
    createChainSet(queries, queries_path);

    sort(queries.begin(), queries.end(), compareChainByLengthAsc);

    std::shared_ptr<Kmers> kmers = createKmers(kmer_length, score_threshold,
        score_matrix);

    uint64_t database_cells = 0;
    std::shared_ptr<Reader> reader = createChainSetPartInitialize(database_path);

    ChainEntrySet entries(queries.size());
    std::vector<MutexPtr> entry_mutexes;
    for (uint32_t i = 0; i < queries.size(); ++i) {
        entry_mutexes.push_back(MutexPtr(new std::mutex()));
    }

    Timer timer;
    while (true) {

        ChainSet database_part;
        auto status = createChainSetPart(database_part, reader, kDatabasePartSize);

        std::vector<uint32_t> tasks;
        preprocDatabase(tasks, database_part, thread_pool->num_threads());

        timer.start();

        std::vector<std::future<void>> thread_futures;

        for (uint32_t i = 0; i < tasks.size() - 1; ++i) {
            thread_futures.emplace_back(thread_pool->submit(scoreChains,
                std::ref(entries), std::ref(entry_mutexes), max_candidates,
                std::ref(queries), std::ref(database_part), tasks[i],
                tasks[i + 1], kmers));
        }

        for (const auto& it: thread_futures) {
            it.wait();
        }

        timer.stop();

        for (const auto& it: database_part) {
            database_cells += it->length();
        }

        if (status == false) {
            break;
        }
    }

    timer.print("database", "search-werk");

    dst.clear();
    dst.resize(queries.size());

    for (uint32_t i = 0; i < queries.size(); ++i) {

        dst[i].reserve(entries[i].size());

        for (const auto& it: entries[i]) {
            dst[i].emplace_back(it.chain_idx());
        }

        std::sort(dst[i].begin(), dst[i].end());
    }

    return database_cells;
}
