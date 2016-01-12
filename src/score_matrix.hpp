/*!
 * @file score_matrix.hpp
 *
 * @brief ScoreMatrix class header file
 */

#pragma once

#include <stdint.h>
#include <memory>
#include <vector>
#include <string>

class ScoreMatrix;

enum class ScoreMatrixType {
    kBlosum45,
    kBlosum50,
    kBlosum62,
    kBlosum80,
    kBlosum90,
    kPam30,
    kPam70,
    kPam250
};

std::unique_ptr<ScoreMatrix> createScoreMatrix(ScoreMatrixType type, int32_t gap_open,
    int32_t gap_extend);

class ScoreMatrix {
public:

    ~ScoreMatrix() = default;

    int* data() {
        return matrix_.data();
    }

    const int* data() const {
        return matrix_.data();
    }

    ScoreMatrixType type() const {
        return type_;
    }

    int32_t gap_open() const {
        return gap_open_;
    }

    int32_t gap_extend() const {
        return gap_extend_;
    }

    std::string scorerName() const;

    int score(uint32_t row, uint32_t column) const;

    friend std::unique_ptr<ScoreMatrix> createScoreMatrix(ScoreMatrixType type,
        int32_t gap_open, int32_t gap_extend);

    static uint32_t num_rows_;
    static uint32_t num_columns_;
    static uint32_t size_;

private:

    ScoreMatrix(ScoreMatrixType type, int32_t gap_open, int32_t gap_extend);
    ScoreMatrix(const ScoreMatrix&) = delete;
    const ScoreMatrix& operator=(const ScoreMatrix&) = delete;

    ScoreMatrixType type_;
    int32_t gap_open_;
    int32_t gap_extend_;
    std::vector<int> matrix_;
};
