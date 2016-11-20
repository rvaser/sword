/*!
 * @file hash.cpp
 *
 * @brief Hash class source file
 */

#include <assert.h>
#include <queue>

#include "chain.hpp"
#include "kmers.hpp"
#include "hash.hpp"

Hit::Hit(uint32_t id, uint32_t position)
        : id_(id), position_(position) {
}

std::unique_ptr<Hash> createHash(const ChainSet& chains, uint32_t start,
    uint32_t length, std::shared_ptr<Kmers> kmers) {

    assert(chains.size());
    assert(start < chains.size() && start + length <= chains.size());
    assert(kmers);

    return std::unique_ptr<Hash>(new Hash(chains, start, length, kmers));
}

Hash::Hash(const ChainSet& chains, uint32_t start, uint32_t length,
    std::shared_ptr<Kmers> kmers)
        : starts_((kmers->mode() == 0 ? kProtNumKmers[kmers->kmer_length()] : kNuclNumKmers[kmers->kmer_length()]), 0) {

    for (uint32_t i = start; i < start + length; ++i) {

        auto kmer_vector = createKmerVector(chains[i], kmers->mode(), kmers->kmer_length());

        for (uint32_t j = 0; j < kmer_vector.size(); ++j) {
            ++starts_[kmer_vector[j] + 1];
            if (kmers->mode() == 0 && kmers->has_permutations()) {
                for (const auto& it: kmers->kmer_substitutions(kmer_vector[j])) {
                    ++starts_[it + 1];
                }
            }
        }
    }

    for (uint32_t i = 1; i < starts_.size() - 1; ++i) {
        starts_[i + 1] += starts_[i];
    }

    hits_.resize(starts_[starts_.size() - 1]);
    std::vector<uint32_t> tmp(starts_.begin(), starts_.end());

    for (uint32_t i = start; i < start + length; ++i) {

        auto kmer_vector = createKmerVector(chains[i], kmers->mode(), kmers->kmer_length());

        for (uint32_t j = 0; j < kmer_vector.size(); ++j) {

            auto hit = Hit(i - start, j);
            hits_[tmp[kmer_vector[j]]++] = hit;
            if (kmers->mode() == 0 && kmers->has_permutations()) {
                for (const auto& it: kmers->kmer_substitutions(kmer_vector[j])) {
                    hits_[tmp[it]++] = hit;
                }
            }
        }
    }
}

void Hash::hits(Iterator& start, Iterator& end, uint32_t key) {
    start = hits_.begin() + starts_[key];
    end = hits_.begin() + starts_[key+1];
}
