/*!
 * @file kmers.cpp
 *
 * @brief Kmers class source file
 */

#include <assert.h>

#include "chain.hpp"
#include "score_matrix.hpp"
#include "kmers.hpp"

constexpr uint32_t kProtMaxValue = 25;
constexpr uint32_t kProtBitLength = 5;

std::vector<uint32_t> kDelMask = { 0, 0, 0, 0x7FFF, 0xFFFFF, 0x1FFFFFF };

std::vector<char> kAminoAcids = {
    /* A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y */
    0, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 15, 16, 17, 18, 19, 21, 22, 24
};

static size_t numKmers(size_t kmer_length) {

    size_t num_kmers = 0;

    for (size_t i = 0; i < kmer_length; ++i) {
        num_kmers += kProtMaxValue << (i * kProtBitLength);
    }

    return num_kmers + 1;
}

static void createKmersRecursive(std::vector<std::string>& kmers,
    std::string& current_kmer, uint32_t level) {

    if (level == 0) {
        kmers.emplace_back(current_kmer);
        return;
    }

    for (const auto& aa: kAminoAcids) {
        current_kmer.push_back(aa);
        createKmersRecursive(kmers, current_kmer, level - 1);
        current_kmer.pop_back();
    }
}

std::vector<uint32_t> createKmerVector(const std::unique_ptr<Chain>& chain,
    uint32_t kmer_length) {

    if (chain->length() < kmer_length) {
        return std::vector<uint32_t>();
    }

    auto data = chain->data();

    std::vector<uint32_t> res(data.size() - kmer_length + 1);
    uint32_t ptr = 0, kmer = 0, del_mask = kDelMask[kmer_length];

    for (uint32_t i = 0; i < kmer_length; ++i) {
        kmer = (kmer << kProtBitLength) | data[i];
    }
    res[ptr++] = kmer;

    for (uint32_t i = kmer_length; i < data.size(); ++i) {
        kmer = ((kmer << kProtBitLength) | data[i]) & del_mask;
        res[ptr++] = kmer;
    }

    return res;
}

std::unique_ptr<Kmers> createKmers(uint32_t kmer_length, uint32_t score_threshold,
    std::shared_ptr<ScoreMatrix> score_matrix) {

    assert(kmer_length > 2);
    assert(kmer_length < 6);
    assert(score_matrix);

    return std::unique_ptr<Kmers>(new Kmers(kmer_length, score_threshold, score_matrix));
}

Kmers::Kmers(uint32_t kmer_length, uint32_t score_threshold, std::shared_ptr<ScoreMatrix> score_matrix) {

    kmer_length_ = kmer_length;
    data_.resize(numKmers(kmer_length));

    if (score_threshold > 0) {
        if (kmer_length_ == 3) {
            createSubstitutionsLong(score_threshold, score_matrix);
        } else {
            createSubstitutionsShort(score_threshold, score_matrix);
        }
    }
}

const std::vector<uint32_t>& Kmers::kmer_substitutions(uint32_t kmer) const {
    return data_[kmer];
}

void Kmers::createSubstitutionsShort(int score_threshold, std::shared_ptr<ScoreMatrix> score_matrix) {

    std::vector<std::string> kmers;
    std::string starting_kmer = "";
    createKmersRecursive(kmers, starting_kmer, kmer_length_);

    for (const auto& kmer_a: kmers) {
        for (size_t i = 0; i < kmer_length_; ++i) {

            auto kmer_b = kmer_a;

            for (const auto& aa: kAminoAcids) {

                if (kmer_a[i] == aa) {
                    continue;
                }

                kmer_b[i] = aa;

                int score = 0;
                for (uint32_t j = 0; j < kmer_length_; ++j) {
                    score += score_matrix->score(kmer_a[j], kmer_b[j]);
                }

                if (score >= score_threshold) {
                    data_[kmer_code(kmer_a)].emplace_back(kmer_code(kmer_b));
                }
            }
        }
    }
}

void Kmers::createSubstitutionsLong(int score_threshold, std::shared_ptr<ScoreMatrix> score_matrix) {

    std::vector<std::string> kmers;
    std::string starting_kmer = "";
    createKmersRecursive(kmers, starting_kmer, kmer_length_);

    for (uint32_t i = 0; i < kmers.size(); ++i) {

        const auto& kmer_a = kmers[i];

        for (uint32_t j = i + 1; j < kmers.size(); ++j) {

            const auto& kmer_b = kmers[j];

            int score = 0;

            for (uint32_t k = 0; k < kmer_length_; ++k) {
                score += score_matrix->score(kmer_a[k], kmer_b[k]);
            }

            if (score >= score_threshold) {
                auto a = kmer_code(kmer_a);
                auto b = kmer_code(kmer_b);
                data_[a].emplace_back(b);
                data_[b].emplace_back(a);
            }
        }
    }
}

uint32_t Kmers::kmer_code(const std::string& kmer) const {

    uint32_t code = 0;

    for (const auto& it: kmer) {
        code <<= kProtBitLength;
        code |= it;
    }

    return code;
}
