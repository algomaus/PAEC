/*
 * BAMIterator.h
 *
 *  Created on: Jan 31, 2017
 *      Author: sarah
 */

#pragma once

#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <string>
#include <vector>

#include "FASTQRead.h"

class FASTQModifiedIterator {
public:
	FASTQModifiedIterator();
	FASTQModifiedIterator(const std::string &readsFileName);
	FASTQRead next();
	std::vector<FASTQRead> next(size_t numReads);
	bool hasReadsLeft();
	size_t numReadsLeft();
	double progress();
private:
	size_t countNumberOfReads(const std::string &readsFileName);
	std::ifstream infileReads;
	size_t readsLeft; size_t numReadsTotal;
};
