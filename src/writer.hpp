/*! copyright DP
 * @file writer.hpp
 *
 * @brief Alignment output header file
 */

#pragma once

#include <cstdio>
#include <vector>
#include <string>
#include <fstream>

class Chain;
class ScoreMatrix;
class Alignment;
class Writer;

using ChainSet = std::vector<std::unique_ptr<Chain>>;
using AlignmentSet = std::vector<std::unique_ptr<Alignment>>;

enum class OutputType {
  kBm0, // BLAST m0
  kBm8, // BLAST m8
  kBm9, // BLAST m9
};

std::unique_ptr<Writer> createWriter(const std::string& path, OutputType format,
  std::shared_ptr<ScoreMatrix> scorer);

class Writer {
 public:
  ~Writer();

  void write_alignments(const AlignmentSet& alignments, const ChainSet& queries,
    const ChainSet& database);

  friend std::unique_ptr<Writer> createWriter(const std::string& path, OutputType format,
    std::shared_ptr<ScoreMatrix> scorer);

 private:
  Writer(FILE* output_file, OutputType format, std::shared_ptr<ScoreMatrix>& scorer);

  Writer(const Writer&) = delete;
  const Writer& operator=(const Writer&) = delete;

  void write_bm0(const AlignmentSet& alignments, const ChainSet& queries,
    const ChainSet& database);

  void write_bm8(const AlignmentSet& alignments, const ChainSet& queries,
    const ChainSet& database);

  void write_bm9(const AlignmentSet& alignments, const ChainSet& queries,
    const ChainSet& database);

  FILE* output_file_;
  OutputType format_;
  std::shared_ptr<ScoreMatrix>& scorer_;
};
