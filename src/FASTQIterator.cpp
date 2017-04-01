/*
 * FASTQIterator.cpp
 *
 *  Created on: Feb 6, 2017
 *      Author: sarah
 */

#include "FASTQIterator.h"

size_t FASTQIterator::countNumberOfReads(const std::string &readsFilename) {
	size_t numReads = 0;
	seqan::SeqFileIn seqFileInCounting;
	if (!open(seqFileInCounting, readsFilename.c_str())) {
		std::cerr << "ERROR: Could not open " << readsFilename << std::endl;
	}
	seqan::CharString id;
	seqan::Dna5String seq;
	seqan::CharString qual;
	while (!atEnd(seqFileInCounting)) {
		seqan::readRecord(id, seq, qual, seqFileInCounting);
		numReads++;
	}

	std::cout << "Finished counting reads. There are " << numReads << " reads in total.\n";
	return numReads;
}

FASTQIterator::FASTQIterator() {
	numReadsTotal = 0;
	readsLeft = 0;
}

FASTQIterator::FASTQIterator(const std::string &readsFilename) {
	numReadsTotal = countNumberOfReads(readsFilename);
	readsLeft = numReadsTotal;

	// Open input file, BamFileIn can read SAM and BAM files.
	if (!open(seqFileIn, readsFilename.c_str())) {
		std::cerr << "ERROR: Could not open " << readsFilename << std::endl;
	}
}

bool FASTQIterator::hasReadsLeft() {
	return (readsLeft > 0);
}

size_t FASTQIterator::numReadsLeft() {
	return readsLeft;
}


double FASTQIterator::progress() {
	return ((numReadsTotal - readsLeft) / ((double) numReadsTotal)) * 100;
}

FASTQRead FASTQIterator::next() {
	FASTQRead fastqRead;
	seqan::CharString id;
	seqan::Dna5String seq;
	seqan::CharString qual;

	if (!atEnd(seqFileIn)) {
		seqan::readRecord(id, seq, qual, seqFileIn);
		fastqRead.id = toCString(id);
		seqan::CharString sequenceAsCharString = seq;
		fastqRead.sequence = toCString(sequenceAsCharString);
		fastqRead.quality = toCString(qual);
		readsLeft--;
	} else {
		throw std::runtime_error("There are no reads left!");
	}
	return fastqRead;
}

std::vector<FASTQRead> FASTQIterator::next(size_t numReads) {
	std::vector<FASTQRead> reads;
	seqan::StringSet<seqan::CharString> ids;
	seqan::StringSet<seqan::Dna5String> seqs;
	seqan::StringSet<seqan::CharString> quals;

	if (!atEnd(seqFileIn)) {
		seqan::readRecords(ids, seqs, quals, seqFileIn, numReads);
		for (size_t i = 0; i < length(ids); ++i) {
			FASTQRead fastqRead;
			fastqRead.id = toCString(ids[i]);
			seqan::CharString sequenceAsCharString = seqs[i];
			fastqRead.sequence = toCString(sequenceAsCharString);
			fastqRead.quality = toCString(quals[i]);
			reads.push_back(fastqRead);
		}
		readsLeft-= reads.size();
	} else {
		throw std::runtime_error("There are no reads left!");
	}
	return reads;
}

