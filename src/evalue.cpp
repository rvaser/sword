/*!
 * @file evalue.cpp
 *
 * @brief EValue class source file
 */

#include <math.h>
#include <vector>

#include "chain.hpp"
#include "reader.hpp"
#include "score_matrix.hpp"
#include "evalue.hpp"

struct EValueConstants {
    int32_t gap_open;
    int32_t gap_extend;
    double lambda;
    double K;
    double H;
    double a;
    double C;
    double alpha;
    double sigma;
};

std::vector<EValueConstants> kEValueConstants = {
    { -1, -1, 0.3176, 0.134, 0.4012, 0.7916, 0.623757, 4.964660, 4.964660 },
    { 11, 2, 0.297, 0.082, 0.27, 1.1, 0.641766, 12.673800, 12.757600 },
    { 10, 2, 0.291, 0.075, 0.23, 1.3, 0.649362, 16.474000, 16.602600 },
    { 9, 2, 0.279, 0.058, 0.19, 1.5, 0.659245, 22.751900, 22.950000 },
    { 8, 2, 0.264, 0.045, 0.15, 1.8, 0.672692, 35.483800, 35.821300 },
    { 7, 2, 0.239, 0.027, 0.10, 2.5, 0.702056, 61.238300, 61.886000 },
    { 6, 2, 0.201, 0.012, 0.061, 3.3, 0.740802, 140.417000, 141.882000 },
    { 13, 1, 0.292, 0.071, 0.23, 1.2, 0.647715, 19.506300, 19.893100 },
    { 12, 1, 0.283, 0.059, 0.19, 1.5, 0.656391, 27.856200, 28.469900 },
    { 11, 1, 0.267, 0.041, 0.14, 1.9, 0.669720, 42.602800, 43.636200 },
    { 10, 1, 0.243, 0.024, 0.10, 2.5, 0.693267, 83.178700, 85.065600 },
    { 9, 1, 0.206, 0.010, 0.052, 4.0, 0.731887, 210.333000, 214.842000 }
};

std::unique_ptr<EValue> createEValue(uint64_t database_cells,
    std::shared_ptr<ScoreMatrix> scorer) {

    return std::unique_ptr<EValue>(new EValue(database_cells, scorer));
}

EValue::EValue(uint64_t database_cells, std::shared_ptr<ScoreMatrix> scorer) {

    uint32_t index_default = 0;
    uint32_t index = 0;

    auto gap_open = scorer->gap_open();
    auto gap_extend = scorer->gap_extend();

    if (scorer->type() == ScoreMatrixType::kBlosum62) {
        for (uint32_t i = 0; i < kEValueConstants.size(); ++i) {
            if (gap_open == kEValueConstants[i].gap_open &&
                    gap_extend == kEValueConstants[i].gap_extend) {
                index = i;
                break;
            }
        }
    }

    double alphaUn = kEValueConstants[index_default].alpha;
    double aUn = kEValueConstants[index_default].a;
    double G = gap_open + gap_extend;

    G_ = G;
    aUn_ = aUn;
    alphaUn_ = alphaUn;
    lambda_ = kEValueConstants[index].lambda;
    K_ = kEValueConstants[index].K;
    logK_ = log(K_);
    H_ = kEValueConstants[index].H;
    a_ = kEValueConstants[index].a;
    C_ = kEValueConstants[index].C;
    alpha_ = kEValueConstants[index].alpha;
    sigma_ = kEValueConstants[index].sigma;
    b_ = 2.0 * G * (aUn_ - a_);
    beta_ = 2.0 * G * (alphaUn_ - alpha_);
    tau_ = 2.0 * G * (alphaUn_ - sigma_);

    length_ = database_cells;

    fprintf(stderr, "Using: lambda = %.3lf, K = %.3lf, H = %.3lf\n",
        lambda_, K_, H_);
}

double EValue::calculate(int32_t score, uint32_t query_length,
    uint32_t target_length) const {

    /* Code taken from SW# (author Matija Korpar) which was take from BLAST */
    int y_ = score;
    int m_ = query_length;
    int n_ = target_length;

    /* The pair-wise e-value must be scaled back to db-wise e-value*/
    double db_scale_factor = (double) length_ / (double) n_;

    double k_         = K_;
    double ai_hat_    = a_;
    double bi_hat_    = b_;
    double alphai_hat_= alpha_;
    double betai_hat_ = beta_;
    double sigma_hat_ = sigma_;
    double tau_hat_   = tau_;

    /* Here we consider symmetric matrix only */
    double aj_hat_    = ai_hat_;
    double bj_hat_    = bi_hat_;
    double alphaj_hat_= alphai_hat_;
    double betaj_hat_ = betai_hat_;

    /* This is 1/sqrt(2.0*PI) */
    static double const_val = 0.39894228040143267793994605993438;
    double m_li_y, vi_y, sqrt_vi_y, m_F, P_m_F;
    double n_lj_y, vj_y, sqrt_vj_y, n_F, P_n_F;
    double c_y, p1, p2;
    double area;

    m_li_y = m_ - (ai_hat_*y_ + bi_hat_);
    vi_y = std::max(2.0*alphai_hat_/lambda_, alphai_hat_*y_+betai_hat_);
    sqrt_vi_y = sqrt(vi_y);
    m_F = m_li_y/sqrt_vi_y;
    P_m_F = 0.5 + 0.5 * erf(m_F);
    p1 = m_li_y * P_m_F + sqrt_vi_y * const_val * exp(-0.5*m_F*m_F);

    n_lj_y = n_ - (aj_hat_*y_ + bj_hat_);
    vj_y = std::max(2.0*alphaj_hat_/lambda_, alphaj_hat_*y_+betaj_hat_);
    sqrt_vj_y = sqrt(vj_y);
    n_F = n_lj_y/sqrt_vj_y;
    P_n_F = 0.5 + 0.5 * erf(n_F);
    p2 = n_lj_y * P_n_F + sqrt_vj_y * const_val * exp(-0.5*n_F*n_F);

    c_y = std::max(2.0*sigma_hat_/lambda_, sigma_hat_*y_+tau_hat_);
    area = p1 * p2 + c_y * P_m_F * P_n_F;

    return area * k_ * exp(-lambda_ * y_) * db_scale_factor;
}
