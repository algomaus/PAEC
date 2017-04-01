/*
 * BAMIterator.cpp
 *
 *  Created on: Jan 31, 2017
 *      Author: sarah
 */

#include "BAMIterator.h"

#include <seqan/bam_io.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <iostream>
#include <stdexcept>

#include "ReadWithAlignments.h"

void BAMIterator::countNumberOfReads(const std::string &alignmentFilename) {
	numReadsTotal = 0;
	numReadsTotalMapped = 0;
	seqan::BamFileIn bamFileInReadCounting;
	if (!open(bamFileInReadCounting, alignmentFilename.c_str())) {
		std::cerr << "ERROR: Could not open " << alignmentFilename << std::endl;
	}
	// read and discard header
	seqan::BamHeader header;
	seqan::readHeader(header, bamFileInReadCounting);
	seqan::BamAlignmentRecord record;
	if (!atEnd(bamFileInReadCounting)) {
		seqan::readRecord(record, bamFileInReadCounting);
		currentReadName = record.qName;
		if (!hasFlagUnmapped(record)) {
			numReadsTotalMapped++;
		}
		numReadsTotal++;
	}
	while (!atEnd(bamFileInReadCounting)) {
		seqan::readRecord(record, bamFileInReadCounting);
		if (record.qName != currentReadName) {
			numReadsTotal++;
			if (!hasFlagUnmapped(record)) {
				numReadsTotalMapped++;
			}
			currentReadName = record.qName;
		}
	}

	std::cout << "Finished counting reads. There are " << numReadsTotal << " reads in total.\n";
}

BAMIterator::BAMIterator() {
	numReadsTotal = 0;
	numReadsTotalMapped = 0;
	readsLeft = 0;
}

BAMIterator::BAMIterator(const std::string &alignmentFilename) {
	countNumberOfReads(alignmentFilename);
	readsLeft = numReadsTotal;

	// Open input file, BamFileIn can read SAM and BAM files.
	if (!open(bamFileIn, alignmentFilename.c_str())) {
		std::cerr << "ERROR: Could not open " << alignmentFilename << std::endl;
	}

	// read and discard header
	seqan::BamHeader header;
	seqan::readHeader(header, bamFileIn);
	seqan::BamAlignmentRecord record;

	if (!atEnd(bamFileIn)) {
		seqan::readRecord(record, bamFileIn);
		currentReadName = record.qName;
		records.push_back(record);
	}
}

bool BAMIterator::hasReadsLeft() {
	return (readsLeft > 0);
}

size_t BAMIterator::numReadsLeft() {
	return readsLeft;
}

double BAMIterator::progress() {
	return ((numReadsTotal - readsLeft) / ((double) numReadsTotal)) * 100;
}

// Assumes that the BAM file is sorted by read name
ReadWithAlignments BAMIterator::next() {
	//std::cout << "next called.\n";
	seqan::BamAlignmentRecord record;

	while (!atEnd(bamFileIn)) {
		seqan::readRecord(record, bamFileIn);
		if (record.qName == currentReadName) {
			records.push_back(record);
		} else {
			ReadWithAlignments alignedRead(records);
			readsLeft--;
			currentReadName = record.qName;
			records.clear();
			records.push_back(record);
			return alignedRead;
		}
	}

	if (records.empty()) {
		throw std::runtime_error("The records are empty!");
	}

	ReadWithAlignments alignedRead(records);
	readsLeft--;
	records.clear();

	assert(readsLeft == 0);
	return alignedRead;
}

std::vector<ReadWithAlignments> BAMIterator::next(size_t numReads) {
	std::vector<ReadWithAlignments> res;
	while (res.size() < numReads && readsLeft > 0) {
		res.push_back(next());
	}
	return res;
}

size_t BAMIterator::getNumReadsTotal() {
	return numReadsTotal;
}

size_t BAMIterator::getNumReadsTotalMapped() {
	return numReadsTotalMapped;
}

