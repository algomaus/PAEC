/*
 * ErrorProfileUnit.cpp
 *
 *  Created on: Jan 31, 2017
 *      Author: sarah
 */

#include "OverallErrorProfile.h"

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "../AlignedInformation/CorrectedReadAligned.h"
#include "../AlignedInformation/CorrectionAligned.h"
#include "../Correction.h"

OverallErrorProfile::OverallErrorProfile() {
	totalCount = 0;
	noncorrectBases = 0;
	deletedBases = 0;
	finalized = false;

	for (ErrorType type : errorTypeIterator()) {
		if (type != ErrorType::CORRECT && type != ErrorType::NODEL) {
			counts[type] = 0;
		}
	}

	std::vector<char> bases = { 'A', 'C', 'G', 'T', 'N' };
	for (size_t i = 0; i < 4; ++i) { // without 'N'
		for (size_t j = 0; j < 5; ++j) { // with 'N'
			substitutionMatrix[std::make_pair(bases[i], bases[j])] = 0;
		}
	}

	overallErrorProbCurrent = 0;
	overallErrorProbNext = 0;
}

double OverallErrorProfile::getOverallErrorRateCurrentBase() {
	if (!finalized) {
		return 1.0 - (double) (totalCount - noncorrectBases) / totalCount;
	} else {
		return overallErrorProbCurrent;
	}
}

double OverallErrorProfile::getOverallErrorRateNextGap() {
	if (!finalized) {
		return 1.0 - (double) (totalCount - deletedBases) / totalCount;
	} else {
		return overallErrorProbNext;
	}
}

void OverallErrorProfile::check(const CorrectedRead &corrRead, double acceptProb) {
	finalized = false;

	totalCount += corrRead.originalRead.sequence.size();
	for (Correction corr : corrRead.corrections) {
		if (corr.type == ErrorType::INSERTION) {
			counts[ErrorType::INSERTION]++;
			noncorrectBases++;
		} else if (corr.type == ErrorType::SUB_FROM_A || corr.type == ErrorType::SUB_FROM_C
				|| corr.type == ErrorType::SUB_FROM_G || corr.type == ErrorType::SUB_FROM_T) {
			substitutionMatrix[std::make_pair(corr.correctedBases[0], corr.originalBases[0])]++;
			noncorrectBases++;
		} else { // corr.type is a deletion, a chimeric break or a multidel
			counts[corr.type]++;
			deletedBases++;
		}
	}
}

// TODO: Fix this code duplication issue.
void OverallErrorProfile::checkAligned(const CorrectedReadAligned &corrRead, double acceptProb) {
	finalized = false;

	totalCount += corrRead.originalRead.sequence.size();
	for (CorrectionAligned ca : corrRead.alignedCorrections) {
		Correction corr = ca.correction;
		if (corr.type == ErrorType::INSERTION) {
			counts[ErrorType::INSERTION]++;
			noncorrectBases++;
		} else if (corr.type == ErrorType::SUB_FROM_A || corr.type == ErrorType::SUB_FROM_C
				|| corr.type == ErrorType::SUB_FROM_G || corr.type == ErrorType::SUB_FROM_T) {
			assert(corr.correctedBases[0] != corr.originalBases[0]);
			substitutionMatrix[std::make_pair(corr.correctedBases[0], corr.originalBases[0])]++;
			noncorrectBases++;
		} else { // corr.type is a deletion, a chimeric break or a multidel
			counts[corr.type]++;
			deletedBases++;
		}
	}
}

std::unordered_map<ErrorType, double> OverallErrorProfile::getErrorProbabilitiesFinalized(const std::string &kmer,
		size_t positionInKmer) {
	assert(finalized);
	std::unordered_map<ErrorType, double> overallProb = counts_finalized;
	overallProb[ErrorType::SUB_FROM_A] = substitutionMatrix_finalized[std::make_pair('A', kmer[positionInKmer])];
	overallProb[ErrorType::SUB_FROM_C] = substitutionMatrix_finalized[std::make_pair('C', kmer[positionInKmer])];
	overallProb[ErrorType::SUB_FROM_G] = substitutionMatrix_finalized[std::make_pair('G', kmer[positionInKmer])];
	overallProb[ErrorType::SUB_FROM_T] = substitutionMatrix_finalized[std::make_pair('T', kmer[positionInKmer])];
	return overallProb;
}

std::unordered_map<ErrorType, double> OverallErrorProfile::getErrorProbabilitiesFinalized(const FASTQRead &read,
		size_t positionInRead) {
	if (read.sequence[positionInRead] == '_') {
			throw std::runtime_error("Invalid k-mer!");
		}
	return getErrorProbabilitiesFinalized(read.sequence, positionInRead);
}

std::unordered_map<ErrorType, double> OverallErrorProfile::getKmerErrorProbabilities(const std::string &kmer,
		size_t positionInKmer) {
	if (kmer.find("_") != std::string::npos) {
		throw std::runtime_error("Invalid k-mer!");
	}

	if (finalized) {
		return getErrorProbabilitiesFinalized(kmer, positionInKmer);
	}
	std::unordered_map<ErrorType, double> overallProb;
	for (auto kv : counts) {
		overallProb[kv.first] = (double) kv.second / totalCount;
	}
	overallProb[ErrorType::SUB_FROM_A] = (double) substitutionMatrix[std::make_pair('A', kmer[positionInKmer])]
			/ totalCount;
	overallProb[ErrorType::SUB_FROM_C] = (double) substitutionMatrix[std::make_pair('C', kmer[positionInKmer])]
			/ totalCount;
	overallProb[ErrorType::SUB_FROM_G] = (double) substitutionMatrix[std::make_pair('G', kmer[positionInKmer])]
			/ totalCount;
	overallProb[ErrorType::SUB_FROM_T] = (double) substitutionMatrix[std::make_pair('T', kmer[positionInKmer])]
			/ totalCount;

	assert(noncorrectBases <= totalCount);
	overallProb[ErrorType::CORRECT] = (double) (totalCount - noncorrectBases) / totalCount;
	assert(deletedBases <= totalCount);
	overallProb[ErrorType::NODEL] = (double) (totalCount - deletedBases) / totalCount;

	for (auto kv : overallProb) {
		overallProb[kv.first] = log(overallProb[kv.first]);
	}

	return overallProb;
}

std::unordered_map<ErrorType, double> OverallErrorProfile::getErrorProbabilities(const FASTQRead &read,
		size_t positionInRead) {
	return getKmerErrorProbabilities(read.sequence, positionInRead);
}

void OverallErrorProfile::storeErrorProfile(const std::string &filepath) {
	std::ofstream outfile(filepath, std::ios::binary);
	cereal::BinaryOutputArchive oarchive(outfile);
	oarchive(*this);
}

void OverallErrorProfile::loadErrorProfile(const std::string &filepath, KmerCounter &counter) {
	std::ifstream infile(filepath, std::ios::binary);
	if (!infile.good()) {
		throw std::runtime_error("The file " + filepath + " does not exist!");
	}
	cereal::BinaryInputArchive iarchive(infile);
	OverallErrorProfile oep;
	iarchive(oep);
	totalCount = oep.totalCount;
	noncorrectBases = oep.noncorrectBases;
	deletedBases = oep.deletedBases;
	counts = oep.counts;
	substitutionMatrix = oep.substitutionMatrix;
	finalized = false;
	finalize();
}

void OverallErrorProfile::plotErrorProfile() {
	std::vector<char> bases = { 'A', 'C', 'G', 'T', 'N' };

	assert(totalCount > 0);

	for (auto kv : counts) {
		if (kv.first != ErrorType::SUB_FROM_A && kv.first != ErrorType::SUB_FROM_C && kv.first != ErrorType::SUB_FROM_G
				&& kv.first != ErrorType::SUB_FROM_T)
			std::cout << "P[" << kv.first << "] = " << log((double) kv.second / totalCount) << "\n";
	}

	for (char invalidBase : bases) {
		if (invalidBase != 'A') {
			std::cout << "P[A <- " << invalidBase << "] = "
					<< log((double) substitutionMatrix[std::make_pair('A', invalidBase)] / totalCount) << "\n";
		}
		if (invalidBase != 'C') {
			std::cout << "P[C <- " << invalidBase << "] = "
					<< log((double) substitutionMatrix[std::make_pair('C', invalidBase)] / totalCount) << "\n";
		}
		if (invalidBase != 'G') {
			std::cout << "P[G <- " << invalidBase << "] = "
					<< log((double) substitutionMatrix[std::make_pair('G', invalidBase)] / totalCount) << "\n";
		}
		if (invalidBase != 'T') {
			std::cout << "P[T <- " << invalidBase << "] = "
					<< log((double) substitutionMatrix[std::make_pair('T', invalidBase)] / totalCount) << "\n";
		}
	}

	std::cout << "P[CORRECT] = " << log((double) (totalCount - noncorrectBases) / totalCount) << "\n";
	std::cout << "P[NODEL] = " << log((double) (totalCount - deletedBases) / totalCount) << "\n";
}

void OverallErrorProfile::reset() {
	totalCount = 0;
	noncorrectBases = 0;
	deletedBases = 0;
	finalized = false;
	for (auto kv : counts) {
		counts[kv.first] = 0;
	}
	for (auto kv : substitutionMatrix) {
		substitutionMatrix[kv.first] = 0;
	}
}

void OverallErrorProfile::finalize() {
	if (finalized) {
		return;
	}
	assert(totalCount > 0);

	for (auto kv : counts) {
		counts_finalized[kv.first] = log((double) kv.second / totalCount);
	}
	std::vector<char> bases = { 'A', 'C', 'G', 'T', 'N' };

	for (char invalidBase : bases) {
		if (invalidBase != 'A') {
			substitutionMatrix_finalized[std::make_pair('A', invalidBase)] = log(
					(double) substitutionMatrix[std::make_pair('A', invalidBase)] / totalCount);
		}
		if (invalidBase != 'C') {
			substitutionMatrix_finalized[std::make_pair('C', invalidBase)] = log(
					(double) substitutionMatrix[std::make_pair('C', invalidBase)] / totalCount);
		}
		if (invalidBase != 'G') {
			substitutionMatrix_finalized[std::make_pair('G', invalidBase)] = log(
					(double) substitutionMatrix[std::make_pair('G', invalidBase)] / totalCount);
		}
		if (invalidBase != 'T') {
			substitutionMatrix_finalized[std::make_pair('T', invalidBase)] = log(
					(double) substitutionMatrix[std::make_pair('T', invalidBase)] / totalCount);
		}
	}

	assert(noncorrectBases <= totalCount);
	counts_finalized[ErrorType::CORRECT] = log((double) (totalCount - noncorrectBases) / totalCount);
	assert(deletedBases <= totalCount);
	counts_finalized[ErrorType::NODEL] = log((double) (totalCount - deletedBases) / totalCount);

	overallErrorProbCurrent = 1.0 - (double) (totalCount - noncorrectBases) / totalCount;
	overallErrorProbNext = 1.0 - (double) (totalCount - deletedBases) / totalCount;

	finalized = true;
}
