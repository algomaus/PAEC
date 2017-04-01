/*
 * BAMIterator.h
 *
 *  Created on: Jan 31, 2017
 *      Author: sarah
 */

#pragma once

#include <seqan/bam_io.h>
#include <seqan/sequence.h>
#include <string>
#include <vector>

#include "ReadWithAlignments.h"

class BAMIterator {
public:
	BAMIterator();
	BAMIterator(const std::string &alignmentFilename);
	ReadWithAlignments next();
	std::vector<ReadWithAlignments> next(size_t numReads);
	bool hasReadsLeft();
	size_t numReadsLeft();
	double progress();
	size_t getNumReadsTotal();
	size_t getNumReadsTotalMapped();
private:
	void countNumberOfReads(const std::string &alignmentFilename);
	seqan::BamFileIn bamFileIn;
	size_t readsLeft; size_t numReadsTotal; size_t numReadsTotalMapped;

	std::vector<seqan::BamAlignmentRecord> records;seqan::CharString currentReadName;
};
