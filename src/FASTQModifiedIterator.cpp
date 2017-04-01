/*
 * FASTQModifiedIterator.cpp
 *
 *  Created on: Feb 6, 2017
 *      Author: sarah
 */

#include "FASTQModifiedIterator.h"
#include <fstream>

size_t FASTQModifiedIterator::countNumberOfReads(const std::string &readsFilename) {
	size_t numReads = 0;
	std::ifstream infile(readsFilename);
	std::string line;
	while (std::getline(infile, line)) {
		if (line.size() > 0) {
			if (line.find_first_not_of("ACGT_") == std::string::npos) {
				numReads++;
			}
		}
	}
	std::cout << "Finished counting reads. There are " << numReads << " reads in total.\n";
	infile.close();
	return numReads;
}

FASTQModifiedIterator::FASTQModifiedIterator() {
	numReadsTotal = 0;
	readsLeft = 0;
}

FASTQModifiedIterator::FASTQModifiedIterator(const std::string &readsFilename) {
	numReadsTotal = countNumberOfReads(readsFilename);
	readsLeft = numReadsTotal;

	infileReads = std::ifstream (readsFilename);

	if (!infileReads.good()) {
		std::cerr << "ERROR: Could not open " << readsFilename << std::endl;
	}
}

bool FASTQModifiedIterator::hasReadsLeft() {
	return (readsLeft > 0);
}

size_t FASTQModifiedIterator::numReadsLeft() {
	return readsLeft;
}


double FASTQModifiedIterator::progress() {
	return ((numReadsTotal - readsLeft) / ((double) numReadsTotal)) * 100;
}

FASTQRead FASTQModifiedIterator::next() {
	if (readsLeft <= 0) {
		throw std::runtime_error("There are no reads left!");
	}
	FASTQRead fastqRead;
	std::string id;
	std::string seq;
	std::string junk;
	std::string qual;

	std::getline(infileReads, id);
	std::getline(infileReads, seq);
	std::getline(infileReads, junk);
	std::getline(infileReads, qual);

	fastqRead.id = id;
	fastqRead.sequence = seq;
	fastqRead.quality = qual;

	readsLeft--;
	return fastqRead;
}

std::vector<FASTQRead> FASTQModifiedIterator::next(size_t numReads) {
	std::vector<FASTQRead> reads;
	for (size_t i = 0; i < numReads; ++i) {
		if (readsLeft == 0) {
			break;
		}
		reads.push_back(next());
	}
	return reads;
}

