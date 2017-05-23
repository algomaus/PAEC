/*
 * ErrorProfileUnit.h
 *
 *  Created on: Jan 31, 2017
 *      Author: sarah
 */

#pragma once

#include <stddef.h>
#include <cassert>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "../external/cereal/archives/binary.hpp"

#include "../AlignedInformation/CorrectedReadAligned.h"
#include "../CorrectedRead.h"
#include "../ErrorType.h"
#include "../FASTQRead.h"

class KmerCounter;

class ErrorProfileUnit {
public:
	virtual ~ErrorProfileUnit() {};

	virtual std::vector<std::unordered_map<ErrorType, double> > getReadErrorProbabilities(const FASTQRead &read) {
		std::vector<std::unordered_map<ErrorType, double> > res;
		for (size_t i = 0; i < read.sequence.size(); ++i) {
			res.push_back(getErrorProbabilities(read, i));
		}
		return res;
	}

	virtual std::vector<std::unordered_map<ErrorType, double> > getKmerErrorProbabilities(const std::string &kmer) {
			std::vector<std::unordered_map<ErrorType, double> > res;
			for (size_t i = 0; i < kmer.size(); ++i) {
				res.push_back(getKmerErrorProbabilities(kmer, i));
			}
			return res;
		}

	virtual std::vector<std::unordered_map<ErrorType, double> > getReadErrorProbabilitiesPartial(const FASTQRead &read, size_t from, size_t to) {
			std::vector<std::unordered_map<ErrorType, double> > res;
			assert(from < read.sequence.size());
			assert(to < read.sequence.size());
			for (size_t i = from; i <= to; ++i) {
				res.push_back(getErrorProbabilities(read, i));
			}
			return res;
		}

	virtual void learnErrorProfileFromFiles(const std::string &correctionsFile, double acceptProb = 1.0) {
		reset();
		std::ifstream infile;
		infile.open(correctionsFile, std::ios::binary);
		if (!infile.good()) {
			throw std::runtime_error("The file " + correctionsFile + " does not exist!");
		}
		size_t minProgress = 1;
		size_t n;
		CorrectedRead cr;
		cereal::BinaryInputArchive iarchive(infile);
		iarchive(n);
		for (size_t i = 0; i < n; ++i) {
			iarchive(cr);
			check(cr, acceptProb);
			double progress = (double) i * 100 / n;
			if (progress >= minProgress) {
				std::cout << progress << "%\n";
				minProgress++;
			}
		}
		infile.close();
		finalize();
	}

	virtual void learnErrorProfileFromFilesAligned(const std::string &correctionsFile, double acceptProb = 1.0) {
		reset();
		std::ifstream infile;
		infile.open(correctionsFile, std::ios::binary);
		if (!infile.good()) {
			throw std::runtime_error("The file " + correctionsFile + " does not exist!");
		}
		size_t minProgress = 0;
		unsigned long long n;
		CorrectedReadAligned cra;
		cereal::BinaryInputArchive iarchive(infile);
		iarchive(n);
		for (unsigned long long i = 0; i < n; ++i) {
			iarchive(cra);
			checkAligned(cra, acceptProb);
			double progress = (double) i * 100 / n;
			if (progress >= minProgress) {
				std::cout << progress << "%\n";
				minProgress++;
			}
		}
		infile.close();
		finalize();
	}

	virtual std::unordered_map<ErrorType, double> getErrorProbabilities(const FASTQRead &read, size_t positionInRead) = 0;
	virtual std::unordered_map<ErrorType, double> getKmerErrorProbabilities(const std::string &kmer, size_t positionInKmer) = 0;
	virtual void loadErrorProfile(const std::string &filepath, KmerCounter &counter) = 0;
	virtual void storeErrorProfile(const std::string &filepath) = 0;
	virtual void plotErrorProfile() = 0;

	virtual void reset() = 0;
	virtual void check(const CorrectedRead &corrRead, double acceptProb = 1.0) = 0;
	virtual void checkAligned(const CorrectedReadAligned &corrRead, double acceptProb = 1.0) = 0;
	virtual void finalize() = 0;
protected:
	virtual std::unordered_map<ErrorType, double> getErrorProbabilitiesFinalized(const FASTQRead &read, size_t positionInRead) = 0;
	virtual std::unordered_map<ErrorType, double> getErrorProbabilitiesFinalized(const std::string &kmer, size_t positionInKmer) = 0;
};
