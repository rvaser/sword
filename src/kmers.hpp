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

constexpr uint32_t kProtMaxValue = 25;
constexpr uint32_t kProtBitLength = 5;
constexpr std::array<uint32_t, 6> kProtNumKmers = { 0, 0, 0, 26427, 845627, 27060027 };
constexpr std::array<uint32_t, 6> kProtDelMask = { 0, 0, 0, 0x7FFF, 0xFFFFF, 0x1FFFFFF };

constexpr uint32_t kNuclMaxValue = 3;
constexpr uint32_t kNuclBitLength = 2;
constexpr std::array<uint32_t, 14> kNuclNumKmers = { 0, 0, 0, 0, 0, 0, 0, 0, 65537, 262145, 1048577, 4194305, 16777217, 67108865 };
constexpr std::array<uint32_t, 14> kNuclDelMask = { 0, 0, 0, 0, 0, 0, 0, 0, 0xFFFF, 0x3FFFF, 0xFFFFF, 0x3FFFFF, 0xFFFFFF, 0x3FFFFFF };

std::unique_ptr<Kmers> createKmers(uint32_t mode, uint32_t kmer_length,
    uint32_t score_threshold, std::shared_ptr<ScoreMatrix> score_matrix);

std::vector<std::pair<uint32_t, uint32_t>> createKmerVector(const std::unique_ptr<Chain>& chain,
    uint32_t mode, uint32_t kmer_length);

uint32_t transformNuclChar(uint32_t c);

class Kmers {
public:

    ~Kmers() = default;

    uint32_t kmer_length() const {
        return kmer_length_;
    }

    uint32_t mode() const {
        return mode_;
    }

    bool has_permutations() const {
        return has_permutations_;
    }

    const std::vector<uint32_t>& kmer_substitutions(uint32_t kmer) const;

    friend std::unique_ptr<Kmers> createKmers(uint32_t mode, uint32_t kmer_length,
        uint32_t score_threshold, std::shared_ptr<ScoreMatrix> score_matrix);

    friend Hash;

private:

    Kmers(uint32_t mode, uint32_t kmer_length, uint32_t score_threshold, std::shared_ptr<ScoreMatrix> score_matrix);
    Kmers(const Kmers&) = delete;
    const Kmers& operator=(const Kmers&) = delete;

    void createSubstitutionsShort(int score_threshold, std::shared_ptr<ScoreMatrix> score_matrix);

    void createSubstitutionsLong(int score_threshold, std::shared_ptr<ScoreMatrix> score_matrix);

    uint32_t prot_kmer_code(const std::string& kmer) const;

    uint32_t mode_; // 0 - protein, 1 - dna
    uint32_t kmer_length_;
    std::vector<std::vector<uint32_t>> data_;
    bool has_permutations_;
};
