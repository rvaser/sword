/*!
 * @file evalue.hpp
 *
 * @brief EValue class header file
 */

#pragma once

#include <stdlib.h>
#include <stdint.h>
#include <memory>

enum class ScoreMatrixType;
class EValue;

std::unique_ptr<EValue> createEValue(uint64_t database_cells,
    std::shared_ptr<ScoreMatrix> scorer);

class EValue {
public:

    ~EValue() = default;

    double calculate(int32_t score, uint32_t query_length, uint32_t target_length) const;

    friend std::unique_ptr<EValue> createEValue(uint64_t database_cells,
        std::shared_ptr<ScoreMatrix> scorer);

private:

    EValue(uint64_t database_cells, std::shared_ptr<ScoreMatrix> scorer);

    EValue(const EValue&) = delete;
	const EValue& operator=(const EValue&) = delete;

    double lambda_;
    double K_;
    double logK_;
    double H_;
    double a_;
    double C_;
    double alpha_;
    double sigma_;
    double b_;
    double beta_;
    double tau_;
    double G_;
    double aUn_;
    double alphaUn_;
    uint64_t length_;
};
