/*!
 * @file kmers.hpp
 *
 * @brief Kmers class header file
 */

#pragma once

#include <string>
#include <memory>
#include <vector>

class ScoreMatrix;
class Chain;
class Hash;
class Kmers;

std::unique_ptr<Kmers> createKmers(uint32_t kmer_length, uint32_t score_threshold,
    std::shared_ptr<ScoreMatrix> score_matrix);

std::vector<uint32_t> createKmerVector(const std::unique_ptr<Chain>& chain,
    uint32_t kmer_length);

class Kmers {
public:

    ~Kmers() = default;

    uint32_t kmer_length() const {
        return kmer_length_;
    }

    const std::vector<uint32_t>& kmer_substitutions(uint32_t kmer) const;

    const std::vector<uint32_t>& kmer_substitutions(const std::string& kmer) const {
        return kmer_substitutions(kmer_code(kmer));
    }

    friend std::unique_ptr<Kmers> createKmers(uint32_t kmer_length, uint32_t score_threshold,
        std::shared_ptr<ScoreMatrix> score_matrix);

    friend Hash;

private:

    Kmers(uint32_t kmer_length, uint32_t score_threshold, std::shared_ptr<ScoreMatrix> score_matrix);
    Kmers(const Kmers&) = delete;
    const Kmers& operator=(const Kmers&) = delete;

    void createSubstitutionsShort(int score_threshold, std::shared_ptr<ScoreMatrix> score_matrix);

    void createSubstitutionsLong(int score_threshold, std::shared_ptr<ScoreMatrix> score_matrix);

    uint32_t kmer_code(const std::string& kmer) const;

    uint32_t kmer_length_;
    std::vector<std::vector<uint32_t>> data_;
};
