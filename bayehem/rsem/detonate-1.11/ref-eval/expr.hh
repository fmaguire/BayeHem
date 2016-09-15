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
#include "fasta.hh"
#include "lazycsv.hh"

typedef std::vector<double> expr;

void read_rsem_expr(expr& tau,
                    const std::string& filename,
                    const fasta& fa)
{
  boost::shared_ptr<std::ifstream> ifs = open_or_throw(filename);
  std::string line;
  lazycsv<8, '\t'> lc;

  try {

    // Check that the header is valid.
    getline(*ifs, line);
    lc.parse_line(line);
    if (lc.at<std::string>(0) != "transcript_id")
      throw std::runtime_error("Invalid header (first column should be transcript_id)");
    if (lc.at<std::string>(5) != "TPM")
      throw std::runtime_error("Invalid header (sixth column should be TPM)");

    // Extract the expression.
    std::vector<bool> seen(fa.card, false);
    while (getline(*ifs, line)) {
      // Parse the line.
      lc.parse_line(line);
      // Check the sequence name and get the corresponding idx.
      std::string name = lc.at<std::string>(0);
      if (fa.names_to_idxs.count(name) == 0)
        throw std::runtime_error("Sequence name " + name + " is not found in the corresponding fasta file.");
      size_t idx = fa.names_to_idxs.find(name)->second;
      // Check that the idx has not already been seen.
      if (seen[idx])
        throw std::runtime_error("Duplicate sequence name " + name);
      // Actually extract the expression.
      tau[idx] = lc.at<double>(5) / 1000000.0;
      seen[idx] = true;
    }

    // Check that all idxs were seen.
    for (size_t i = 0; i < fa.card; ++i)
      if (!seen[i])
        throw std::runtime_error("No expression for sequence name " + fa.names[i]);

  } catch (const std::runtime_error& x) {
    throw std::runtime_error("Can't parse " + filename + ": " + x.what());
  }
}
