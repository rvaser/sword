/*! copyright DP
 * @file writer.hpp
 *
 * @brief Alignment output source file
 */

#include <cstdio>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>

#include "chain.hpp"
#include "score_matrix.hpp"
#include "database_alignment.hpp"
#include "writer.hpp"

constexpr int kMatch = 0;  // a match
constexpr int kDeletion = 1;  // deletion from query (insertion to target)
constexpr int kInsertion = 2;  // insertion to query (deletion from target)
constexpr int kMismatch = 3;  // mismatch

std::unique_ptr<Writer> createWriter(const std::string& path, OutputType format,
  std::shared_ptr<ScoreMatrix> scorer) {

  auto output_file = path.empty() ? stdout : fopen(path.c_str(), "w");

  return std::unique_ptr<Writer>(new Writer(output_file, format, scorer));
}

Writer::Writer(FILE* output_file, OutputType format, std::shared_ptr<ScoreMatrix>& scorer)
    : output_file_(output_file), format_(format), scorer_(scorer) {
}

Writer::~Writer() {
  if (output_file_ != stdout) {
    fclose(output_file_);
  }
}

void Writer::write_alignments(const AlignmentSet& alignments, const ChainSet& queries,
  const ChainSet& database) {

    switch (format_) {
      case OutputType::kBm0:
        write_bm0(alignments, queries, database);
        break;

      case OutputType::kBm8:
        write_bm8(alignments, queries, database);
        break;

      case OutputType::kBm9:
        write_bm9(alignments, queries, database);
        break;

      default:
        fprintf(stderr, "Non-existent output type!\n");
        break;
    }
}

void Writer::write_bm0(const AlignmentSet& alignments, const ChainSet& queries,
  const ChainSet& database) {

  if (alignments.size() == 0) {
    fprintf(output_file_, "No alignments found\n");
    return;
  }

  auto query_id = alignments[0]->query_id();
  auto& query = queries[query_id];
  const auto& query_seq = query->data();

  fprintf(output_file_, "Query= %s\n", query->name().c_str());
  fprintf(output_file_, "Length=%zu\n\n", query->length());
  fprintf(output_file_, "Sequences producing significant alignments:");
  fprintf(output_file_, "%27.27s", "Score");
  fprintf(output_file_, "%10.10s\n\n", "Evalue");

  for (const auto& alignment : alignments) {
    auto target_id = alignment->target_id();
    auto& target = database[target_id];

    auto& name = target->name();
    auto score = alignment->score();
    auto eval = alignment->evalue();

    if (name.size() > 64) {
      fprintf(output_file_, "     %.64s...%10d%10.0e\n", name.c_str(), score, eval);
    } else {
      fprintf(output_file_, "     %.67s%10d%10.0e\n", name.c_str(), score, eval);
    }
  }

  fprintf(output_file_, "\n");

  for (const auto& alignment : alignments) {
    auto target_id = alignment->target_id();
    auto& target = database[target_id];
    auto& target_name = target->name();
    const auto& target_seq = target->data();
    auto align_target_start = alignment->target_begin();
    auto align_query_start = alignment->query_begin();

    fprintf(output_file_, ">%s\n", target_name.c_str());
    fprintf(output_file_, "Length=%zu\n\n", target->length());
    fprintf(output_file_, " Score = %d,", alignment->score());
    fprintf(output_file_, " Expect = %.0e\n", alignment->evalue());

    const auto& alignment_str = alignment->alignment();
    int identities = 0;
    int positives = 0;
    int gaps = 0;
    auto alignment_len = alignment_str.size();

    auto query_ptr = align_query_start;
    auto target_ptr = align_target_start;

    for (uint32_t j = 0; j < alignment_len; ++j) {
      int align_res = alignment_str[j];

      if (align_res == kMatch) {
        identities++;
        positives++;
        query_ptr++;
        target_ptr++;
      } else if (align_res == kDeletion) {
        gaps++;
        query_ptr++;
      } else if (align_res == kInsertion) {
        gaps++;
        target_ptr++;
      } else {
        if (scorer_->score(query_seq[query_ptr], target_seq[target_ptr]) > 0) {
          positives++;
        }
        query_ptr++;
        target_ptr++;
      }
    }

    int idn_pct = static_cast<int>(floor(identities*100.f/alignment_len));
    int pos_pct = static_cast<int>(floor(positives*100.f/alignment_len));
    int gap_pct = static_cast<int>(floor(gaps*100.f/alignment_len));

    fprintf(output_file_, " Identities = %d/%lu (%d%%),", identities, alignment_len, idn_pct);
    fprintf(output_file_, " Positives = %d/%lu (%d%%),", positives, alignment_len, pos_pct);
    fprintf(output_file_, " Gaps = %d/%lu (%d%%)", gaps, alignment_len, gap_pct);
    fprintf(output_file_, "\n\n");

    char query_str[61] = { 0 };
    char target_str[61] = { 0 };
    char markup_str[61] = { 0 };

    auto query_start = align_query_start;
    auto query_end = query_start;
    auto target_start = align_target_start;
    auto target_end = target_start;

    for (uint32_t j = 0; j < alignment_len; ++j) {
      int align_res = alignment_str[j];
      int index = j % 60;

      if (align_res == kMatch) {
        markup_str[index] = query_seq[query_end] + 'A';
        query_str[index] = query_seq[query_end++] + 'A';
        target_str[index] = target_seq[target_end++] + 'A';
      } else if (align_res == kMismatch) {
        if (scorer_->score(query_seq[query_end], target_seq[target_end]) > 0) {
          markup_str[index] = '+';
        } else {
          markup_str[index] = ' ';
        }
        query_str[index] = query_seq[query_end++] + 'A';
        target_str[index] = target_seq[target_end++] + 'A';
      } else if (align_res == kDeletion) {
        markup_str[index] = ' ';
        query_str[index] = query_seq[query_end++] + 'A';
        target_str[index] = '-';
      } else {
        markup_str[index] = ' ';
        query_str[index] = '-';
        target_str[index] = target_seq[target_end++] + 'A';
      }

      if ((j+1) % 60 == 0 || j == alignment_len-1) {
        fprintf(output_file_, "Query  %-6d%s  %d\n"
                              "%12s %s %s\n"
                              "Sbjct  %-6d%s  %d\n\n",
                                query_start+1, query_str, query_end,
                                "", markup_str, "",
                                target_start+1, target_str, target_end);

        query_start = query_end;
        target_start = target_end;

        memset(query_str, 0, 61 * sizeof(char));
        memset(target_str, 0, 61 * sizeof(char));
        memset(markup_str, 0, 61 * sizeof(char));
      }
    }
    fprintf(output_file_, "\n");
  }

  fprintf(output_file_, "Matrix: %s\n", scorer_->scorerName().c_str());
  fprintf(output_file_, "Gap Penalties: Existence: %d, Extension: %d\n",
      scorer_->gap_open(), scorer_->gap_extend());
}

void Writer::write_bm8(const AlignmentSet& alignments, const ChainSet& queries,
  const ChainSet& database) {

  if (alignments.size() == 0) {
    return;
  }

  auto query_id = alignments[0]->query_id();
  auto query_name = queries[query_id]->name();

  auto space_pos = query_name.find(" ", 0);
  if (space_pos != std::string::npos) {
    query_name = query_name.substr(0, space_pos);
  }

  for (const auto& alignment : alignments) {
    int mismatches = 0;
    int gap_openings = 0;
    int matches = 0;
    bool gap = false;

    const auto& alignment_str = alignment->alignment();
    const auto alignment_len = alignment_str.size();

    for (uint32_t j = 0; j < alignment_len; ++j) {
      int alignment_res = alignment_str[j];
      // printf("Alignment res %d\n", alignment_res);
      // return;

      switch (alignment_res) {
        case kMismatch:
          mismatches++;
          gap = false;
          break;
        case kInsertion:
          if (!gap) {
            gap = true;
            gap_openings++;
          }
          break;
        case kDeletion:
          if (!gap) {
            gap = true;
            gap_openings++;
          }
          break;
        default:
          matches++;
          gap = false;
      }
    }

    auto target_id = alignment->target_id();
    auto target_name = database[target_id]->name();

    space_pos = target_name.find(" ", 0);
    if (space_pos != std::string::npos) {
      target_name = target_name.substr(0, space_pos);
    }

    double perc_id = (100.f * matches) / alignment_len;

    fprintf(output_file_, "%s\t%s\t", query_name.c_str(), target_name.c_str());
    fprintf(output_file_, "%.lf\t%lu\t%d\t", perc_id, alignment_len, mismatches);
    fprintf(output_file_, "%d\t%d\t", gap_openings, alignment->query_begin()+1);
    fprintf(output_file_, "%d\t%d\t", alignment->query_end()+1, alignment->target_begin()+1);
    fprintf(output_file_, "%d\t", alignment->target_end()+1);

    auto eval = alignment->evalue();
    auto score = alignment->score();

    if (eval >= 1e-2 && eval < 100) {
      fprintf(output_file_, "%.2lf\t", eval);
    } else {
      fprintf(output_file_, "%.2e\t", eval);
    }

    fprintf(output_file_, "%d", score);
    fprintf(output_file_, "\n");
  }
}

void Writer::write_bm9(const AlignmentSet& alignments, const ChainSet& queries,
  const ChainSet& database) {

  fprintf(output_file_, "# Fields:\n");
  fprintf(output_file_, "Query id,Subject id,%% identity,alignment length,mismatches,"
    "gap openings,q. start,q. end,s. start,s. end,e-value,score\n");

  write_bm8(alignments, queries, database);
}
