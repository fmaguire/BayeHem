// Copyright (c) 2013
// Nathanael Fillmore (University of Wisconsin-Madison)
// nathanae@cs.wisc.edu
//
// This file is part of REF-EVAL.
//
// REF-EVAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// REF-EVAL is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with REF-EVAL.  If not, see <http://www.gnu.org/licenses/>.

#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <string>
#include <vector>
#include <boost/foreach.hpp>
#include <sparsehash/sparse_hash_map>
#include <sparsehash/dense_hash_map>
#include "expr.hh"
#include "fasta.hh"
#include "opts.hh"
#include "skip_Ns.hh"
#include "util.hh"
#include "kmer_key.hh"

namespace re {
namespace kc {

template<typename Number>
struct kmer_info
{
  bool is_present_in_A;
  Number weight_in_B;
  kmer_info()
  : is_present_in_A(false),
    weight_in_B(0.0)
  {}
};

typedef google::sparse_hash_map<kmer_key, kmer_info<double>, kmer_key_hash, kmer_key_equal_to> sparse_double_kmer_map;
typedef google::sparse_hash_map<kmer_key, kmer_info<float>, kmer_key_hash, kmer_key_equal_to> sparse_float_kmer_map;
typedef google::dense_hash_map<kmer_key, kmer_info<double>, kmer_key_hash, kmer_key_equal_to> dense_double_kmer_map;
typedef google::dense_hash_map<kmer_key, kmer_info<float>, kmer_key_hash, kmer_key_equal_to> dense_float_kmer_map;


template<typename Ht>
struct empty_key_initializer
{
  empty_key_initializer(Ht&, size_t)
  {}
};

template<>
struct empty_key_initializer<dense_double_kmer_map>
{
  std::string empty_string;
  const char *empty_key;
  empty_key_initializer(dense_double_kmer_map& ht, size_t kmerlen)
  : empty_string(kmerlen, ' '),
    empty_key(empty_string.c_str())
  {
    ht.set_empty_key(empty_key);
  }
};

template<>
struct empty_key_initializer<dense_float_kmer_map>
{
  std::string empty_string;
  const char *empty_key;
  empty_key_initializer(dense_float_kmer_map& ht, size_t kmerlen)
  : empty_string(kmerlen, ' '),
    empty_key(empty_string.c_str())
  {
    ht.set_empty_key(empty_key);
  }
};

template<typename Ht>
void count_kmers_in_A(
    Ht& ht,
    const std::vector<std::string>& A,
    const std::vector<std::string>& A_rc,
    size_t kmerlen,
    bool strand_specific)
{
  // For each contig a in A:
  //   For each kmer r in a or reverse_complement(a):
  //     Mark r as being present in A.
  size_t num_strands = strand_specific ? 1 : 2;
  for (size_t i = 0; i < A.size(); ++i) {
    //std::cerr << i << " of " << A.size() << "(" << 100.0*i/A.size()
    //          << " percent)" << std::endl;
    for (size_t which = 0; which < num_strands; ++which) {
      const std::string& a = which == 0 ? A[i] : A_rc[i];
      if (a.size() >= kmerlen) {
        const char *beg = a.c_str();
        const char *a_end = a.c_str() + a.size() + 1 - kmerlen;
        beg = skip_Ns(beg, a_end, kmerlen, true);
        for (; beg != a_end; ++beg) {
          beg = skip_Ns(beg, a_end, kmerlen, false);
          if (beg == a_end)
            break;
          ht[beg].is_present_in_A = true;
        }
      }
    }
  }
}

template<typename Ht>
void count_kmers_in_B(
    Ht& ht,
    const std::vector<std::string>& B,
    const std::vector<std::string>& B_rc,
    const std::vector<double>& tau_B,
    size_t kmerlen,
    bool strand_specific)
{
  // For each contig b in B:
  //   For each kmer r in b or reverse_complement(b):
  //     Add weight(b) to weight_in_B(r).
  size_t num_strands = strand_specific ? 1 : 2;
  for (size_t i = 0; i < B.size(); ++i) {
    //std::cerr << i << " of " << B.size() << "(" << 100.0*i/B.size()
    //          << " percent)" << std::endl;
    for (size_t which = 0; which < num_strands; ++which) {
      const std::string& b = which == 0 ? B[i] : B_rc[i];
      if (b.size() >= kmerlen) {
        double c = tau_B[i];
        const char *beg = b.c_str();
        const char *b_end = b.c_str() + b.size() + 1 - kmerlen;
        beg = skip_Ns(beg, b_end, kmerlen, true);
        for (; beg != b_end; ++beg) {
          beg = skip_Ns(beg, b_end, kmerlen, false);
          if (beg == b_end)
            break;
          ht[beg].weight_in_B += c; // relies on default init to 0
        }
      }
    }
  }
}

size_t estimate_hashtable_size(
    const std::vector<std::string>& A,
    const std::vector<std::string>& B,
    size_t kmerlen,
    double hash_table_fudge_factor)
{
  size_t max_entries = 0;
  BOOST_FOREACH(const std::string& a, A)
    if (a.size() >= kmerlen)
      max_entries += static_cast<size_t>(0.5 +
        2 * (a.size() + 1 - kmerlen) / hash_table_fudge_factor);
  BOOST_FOREACH(const std::string& b, B)
    if (b.size() >= kmerlen)
      max_entries += static_cast<size_t>(0.5 +
        2 * (b.size() + 1 - kmerlen) / hash_table_fudge_factor);
  return max_entries;
}

template<typename Ht>
double compute_kmer_recall(const Ht& ht)
{
  typedef typename Ht::value_type X;
  double numer = 0.0, denom = 0.0;
  BOOST_FOREACH(const X& x, ht) {
    const typename Ht::mapped_type& i = x.second;
    if (i.is_present_in_A > 0)
      numer += i.weight_in_B;
    denom += i.weight_in_B;
  }
  return numer / denom;
}

double compute_inverse_compression_rate(
    const opts& o,
    const fasta& A)
{
  size_t num_bases_in_A = std::accumulate(A.lengths.begin(), A.lengths.end(), 0);
  return 1.0 * num_bases_in_A / (o.num_reads * o.readlen);
}

template<typename Ht>
void main_1(
    const opts& o,
    const fasta& A,
    const fasta& B,
    const expr& tau_B)
{
  std::cerr << "Reverse complementing the sequences..." << std::endl;
  std::vector<std::string> A_rc, B_rc;
  transform(A.seqs.begin(), A.seqs.end(), back_inserter(A_rc), reverse_complement);
  transform(B.seqs.begin(), B.seqs.end(), back_inserter(B_rc), reverse_complement);

  size_t max_entries = estimate_hashtable_size(A.seqs, B.seqs, o.kmerlen, o.hash_table_fudge_factor);
  std::cerr << "Initializing the hash table with space for " << max_entries << " entries..." << std::endl;
  Ht ht(max_entries, kmer_key_hash(o.kmerlen), kmer_key_equal_to(o.kmerlen));
  empty_key_initializer<Ht> eki(ht, o.kmerlen);

  std::cerr << "Populating the hash table..." << std::flush;
  count_kmers_in_A(ht, A.seqs, A_rc,        o.readlen, o.strand_specific);
  count_kmers_in_B(ht, B.seqs, B_rc, tau_B, o.readlen, o.strand_specific);
  std::cerr << "done; hash table contains " << ht.size() << " entries." << std::endl;

  std::cerr << "Computing kmer recall, inverse compression rate, and kmer compression scores..." << std::endl;
  double wkr = compute_kmer_recall(ht);
  double icr = compute_inverse_compression_rate(o, A);

  std::cout << "weighted_kmer_recall\t" << wkr << std::endl;
  std::cout << "inverse_compression_rate\t" << icr << std::endl;
  std::cout << "kmer_compression_score\t" << wkr - icr << std::endl;
}

void main(
    const opts& o,
    const fasta& A,
    const fasta& B,
    const expr& tau_B)
{
  if (o.kc || o.paper) {
    if (o.hash_table_type == "sparse" && o.hash_table_numeric_type == "double")
      main_1<sparse_double_kmer_map>(o, A, B, tau_B);
    else if (o.hash_table_type == "sparse" && o.hash_table_numeric_type == "float")
      main_1<sparse_float_kmer_map>(o, A, B, tau_B);
    else if (o.hash_table_type == "dense" && o.hash_table_numeric_type == "double")
      main_1<dense_double_kmer_map>(o, A, B, tau_B);
    else if (o.hash_table_type == "dense" && o.hash_table_numeric_type == "float")
      main_1<dense_float_kmer_map>(o, A, B, tau_B);
    else
      throw std::runtime_error("Unknown hash map type.");
  }
}

} // namespace kc
} // namespace re
