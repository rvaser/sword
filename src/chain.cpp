/*!
 * @file chain.cpp
 *
 * @brief Chain class source file
 */

#include <assert.h>

#include "reader.hpp"
#include "chain.hpp"

std::vector<char> kCoder = {
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1,   0,   1,   2,   3,   4,
     5,   6,   7,   8,   9,  10,  11,  12,  13,  14,
    15,  16,  17,  18,  19,  20,  21,  22,  23,  24,
    25,  -1,  -1,  -1,  -1,  -1,  -1,   0,   1,   2,
     3,   4,   5,   6,   7,   8,   9,  10,  11,  12,
    13,  14,  15,  16,  17,  18,  19,  20,  21,  22,
    23,  24,  25,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1
};

std::unique_ptr<Chain> createChain(uint32_t id, char* name, uint32_t name_length,
    char* data, uint32_t data_length) {

    assert(name_length);
    assert(data_length);

    // remove trailing white spaces
    while (isspace(name[name_length - 1])) {
        --name_length;
    }
    assert(name_length && "name cannot be empty");

    std::string data_;
    data_.reserve(data_length);
    uint32_t data_ptr = 0, valid_data_length = 0;

    for (uint32_t i = 0; i < data_length; ++i) {
        auto c = kCoder[data[i]];
        if (c != -1) {
            data[data_ptr++] = c;
            ++valid_data_length;
        }
    }
    assert(valid_data_length && "no valid chars found");

    return std::unique_ptr<Chain>(new Chain(id, std::string(name, name_length),
        std::string(data, valid_data_length)));
}

void createChainSet(ChainSet& dst, const std::string& path) {

    auto reader = createChainSetPartInitialize(path);
    createChainSetPart(dst, std::move(reader), 0);
}

std::unique_ptr<Reader> createChainSetPartInitialize(const std::string& path) {

    /* maybe in future: check if the chains are cached */
    return createReader(path);
}

bool createChainSetPart(ChainSet& dst, std::shared_ptr<Reader> reader, size_t max_bytes) {

    assert(reader);
    return reader->read_chains(dst, max_bytes);
}

Chain::Chain(uint32_t id, std::string&& name, std::string&& data)
        : id_(id), name_(name), data_(data), valid_regions_(), is_dna_(false) {

    valid_regions_.emplace_back(0, data.size());
}

void Chain::change_protein_to_dna_codes() {

    if (is_dna_ == true) {
        return;
    }
    is_dna_ = true;

    valid_regions_.clear();

    bool is_n_region = false;
    uint32_t begin = 0;
    if (data_[0] == 13) {
        is_n_region = true;
    }

    for (uint32_t i = 0; i < data_.size(); ++i) {
        if (data_[i] != 0 && data_[i] != 2 && data_[i] != 6 && data_[i] != 19 && data_[i] != 20) {
            data_[i] = 13; // replace everything ambiguous to N for database search
            if (is_n_region == false) {
                is_n_region = true;
                valid_regions_.emplace_back(begin, i);
            }
        } else {
            switch (data_[i]) {
                case 0: // A
                    break;
                case 2: // C
                    data_[i] = 1;
                    break;
                case 6: // G
                    data_[i] = 2;
                    break;
                case 19: // T
                    data_[i] = 3;
                    break;
                case 20: // U
                    data_[i] = 3;
                    break;
                default:
                    assert(false && "not valid dna");
                    break;
            }
            if (is_n_region == true) {
                is_n_region = false;
                begin = i;
            }
        }
    }

    if (data_.back() != 13) {
        valid_regions_.emplace_back(begin, data_.size());
    }
}
