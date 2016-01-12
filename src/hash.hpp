/*!
 * @file hash.hpp
 *
 * @brief Hash class header file
 */

#pragma once

#include <stdint.h>
#include <memory>
#include <vector>

class Chain;
class Kmers;
class Hash;

using ChainSet = std::vector<std::unique_ptr<Chain>>;

class Hit {
public:

    Hit() {};
    Hit(uint32_t id, uint32_t position);
    ~Hit() = default;

    uint32_t id() const {
        return id_;
    }

    uint32_t position() const {
        return position_;
    }

private:

    uint32_t id_;
    uint32_t position_;
};

std::unique_ptr<Hash> createHash(const ChainSet& chains, uint32_t start,
    uint32_t length, std::shared_ptr<Kmers> kmers);

class Hash {
public:

    ~Hash() {};

    using Iterator = std::vector<Hit>::iterator;
    void hits(Iterator& start, Iterator& end, uint32_t key);

    friend std::unique_ptr<Hash> createHash(const ChainSet& chains,
        uint32_t start, uint32_t length, std::shared_ptr<Kmers> kmers);

private:

    Hash(const ChainSet& chains, uint32_t start, uint32_t length,
        std::shared_ptr<Kmers> kmers);

    Hash(const Hash&) = delete;
    const Hash& operator=(const Hash&) = delete;

    std::vector<size_t> starts_;
    std::vector<Hit> hits_;
};
