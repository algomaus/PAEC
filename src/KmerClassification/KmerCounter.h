/*
 * KmerCounter.h
 *
 *  Created on: Jan 31, 2017
 *      Author: sarah
 */

#pragma once

#include <sdsl/csa_wt.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/rank_support_v5.hpp>
#include <sdsl/select_support_scan.hpp>
#include <sdsl/wt_huff.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <stddef.h>
#include <memory>
#include <string>
#include <unordered_map>
#include "../ErrorProfile/ErrorProfileUnit.hpp"

#include "../ErrorType.h"

using namespace sdsl;

typedef csa_wt<wt_huff<bit_vector, rank_support_v5<>, select_support_scan<>, select_support_scan<0>>, 1 << 20, 1 << 20> FMIndex;

class KmerCounter {
public:
	KmerCounter(const std::string &filepath);
	size_t countKmer(const std::string &kmer);
	double countKmerApproximate(const std::string &kmer, const std::shared_ptr<ErrorProfileUnit> &errorProfile);
	size_t countKmerNoRC(const std::string &kmer);
	double countKmerNoRCApproximate(const std::string &kmer, const std::shared_ptr<ErrorProfileUnit> &errorProfile);
	void clearBuffers();
	void changeFile(const std::string &filepath);
private:
	std::string reverseComplementString(const std::string &sequence);
	std::string kmerAfterError(const std::string &kmer, ErrorType error, size_t posOfError);
	FMIndex fm_index;

	std::unordered_map<std::string, size_t> buffer;
	std::unordered_map<std::string, size_t> bufferApprox;
};

