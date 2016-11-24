/*!
 * @file kmers.cpp
 *
 * @brief Kmers class source file
 */

#include <assert.h>

#include "chain.hpp"
#include "score_matrix.hpp"
#include "kmers.hpp"

std::vector<char> kAminoAcids = {
    /* A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y */
    0, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 15, 16, 17, 18, 19, 21, 22, 24
};

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

std::vector<std::pair<uint32_t, uint32_t>> createKmerVector(const std::unique_ptr<Chain>& chain,
    uint32_t mode, uint32_t kmer_length) {

    if (chain->length() < kmer_length) {
        return std::vector<std::pair<uint32_t, uint32_t>>();
    }

    auto data = chain->data();

    std::vector<std::pair<uint32_t, uint32_t>> dst;
    uint32_t del_mask = mode == 0 ? kProtDelMask[kmer_length] : kNuclDelMask[kmer_length];
    uint32_t bit_length = mode == 0 ? kProtBitLength : kNuclBitLength;

    for (const auto& it: chain->valid_regions()) {
        uint32_t kmer = 0;
        for (uint32_t i = it.first; i < it.second; ++i) {
            kmer = ((kmer << bit_length) | data[i]) & del_mask;
            if (i - it.first >= kmer_length - 1) {
                dst.emplace_back(i - kmer_length + 1, kmer);
            }
        }
    }

    return dst;
}

std::unique_ptr<Kmers> createKmers(uint32_t mode, uint32_t kmer_length,
    uint32_t score_threshold, std::shared_ptr<ScoreMatrix> score_matrix) {

    assert((mode == 0 && kmer_length > 2 && kmer_length < 6) ||
        (mode == 1 && kmer_length > 7 && kmer_length < 14));
    assert(score_matrix);

    return std::unique_ptr<Kmers>(new Kmers(mode, kmer_length, score_threshold, score_matrix));
}

Kmers::Kmers(uint32_t mode, uint32_t kmer_length, uint32_t score_threshold,
    std::shared_ptr<ScoreMatrix> score_matrix)
        : mode_(mode), kmer_length_(kmer_length), data_(), has_permutations_(false){

    if (mode_ == 0 && score_threshold > 0) {
        has_permutations_ = true;
        data_.resize(kProtNumKmers[kmer_length]);
        if (kmer_length_ == 3) {
            createSubstitutionsLong(score_threshold, score_matrix);
        } else {
            createSubstitutionsShort(score_threshold, score_matrix);
        }
    }
}

const std::vector<uint32_t>& Kmers::kmer_substitutions(uint32_t kmer) const {
    assert(mode_ == 0 && "Nucleotide permutations are disabled!");
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
                    data_[prot_kmer_code(kmer_a)].emplace_back(prot_kmer_code(kmer_b));
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
                auto a = prot_kmer_code(kmer_a);
                auto b = prot_kmer_code(kmer_b);
                data_[a].emplace_back(b);
                data_[b].emplace_back(a);
            }
        }
    }
}

uint32_t Kmers::prot_kmer_code(const std::string& kmer) const {

    uint32_t code = 0;

    for (const auto& it: kmer) {
        code <<= kProtBitLength;
        code |= it;
    }

    return code;
}
