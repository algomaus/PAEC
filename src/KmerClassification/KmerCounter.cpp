/*
 * KmerCounter.cpp
 *
 *  Created on: Jan 31, 2017
 *      Author: sarah
 */

#include "KmerCounter.h"

#include <stdexcept>
#include <fstream>

void KmerCounter::changeFile(const std::string &filepath) {
	clearBuffers();
	std::ifstream infile(filepath);
	if (!infile.good()) {
		throw std::runtime_error("This file does not exist! " + filepath);
	}

	std::string index_suffix = ".fm9";
	std::string index_file = filepath + index_suffix;
	if (!load_from_file(fm_index, index_file)) {
		std::cout << "No FM index found. Constructing FM index..." << std::endl;
		construct(fm_index, filepath, 1); // generate index
		store_to_file(fm_index, index_file); // save it
		std::cout << "Index construction complete, index requires " << size_in_mega_bytes(fm_index) << " MiB."
				<< std::endl;
	}
}

KmerCounter::KmerCounter(const std::string &filepath) {
	std::ifstream infile(filepath);
	if (!infile.good()) {
		throw std::runtime_error("This file does not exist! " + filepath);
	}

	std::string index_suffix = ".fm9";
	std::string index_file = filepath + index_suffix;
	if (!load_from_file(fm_index, index_file)) {
		std::cout << "No FM index found. Constructing FM index..." << std::endl;
		construct(fm_index, filepath, 1); // generate index
		store_to_file(fm_index, index_file); // save it
		std::cout << "Index construction complete, index requires " << size_in_mega_bytes(fm_index) << " MiB."
				<< std::endl;
	}
}

void KmerCounter::clearBuffers() {
	buffer.clear();
	bufferApprox.clear();
}

size_t KmerCounter::countKmer(const std::string &kmer) {
	size_t countOriginal = countKmerNoRC(kmer);
	std::string kmerRC = reverseComplementString(kmer);
	size_t countRC = countKmerNoRC(kmerRC);
	return countOriginal + countRC;
}

double KmerCounter::countKmerApproximate(const std::string &kmer,
		const std::shared_ptr<ErrorProfileUnit> &errorProfile) {
	double countOriginal = countKmerNoRCApproximate(kmer, errorProfile);
	std::string kmerRC = reverseComplementString(kmer);
	double countRC = countKmerNoRCApproximate(kmerRC, errorProfile);
	return countOriginal + countRC;
}

size_t KmerCounter::countKmerNoRC(const std::string &kmer) {
	/*if ((kmer.size() < 17) && (buffer.find(kmer) != buffer.end())) {
		return buffer[kmer];
	}*/
	size_t countOriginal = sdsl::count(fm_index, kmer.begin(), kmer.end());
	/*if (kmer.size() < 17 && countOriginal >= 50) {
		buffer[kmer] = countOriginal;
	}*/
	return countOriginal;
}

std::string KmerCounter::kmerAfterError(const std::string &kmer, ErrorType error, size_t posOfError) {
	std::string res = kmer;
	if (error == ErrorType::SUB_FROM_A) {
		res[posOfError] = 'A';
	} else if (error == ErrorType::SUB_FROM_C) {
		res[posOfError] = 'C';
	} else if (error == ErrorType::SUB_FROM_G) {
		res[posOfError] = 'G';
	} else if (error == ErrorType::SUB_FROM_T) {
		res[posOfError] = 'T';
	} else if (error == ErrorType::INSERTION) {
		res = kmer.substr(0, posOfError) + kmer.substr(posOfError + 1, kmer.size());
	} else if (error == ErrorType::DEL_OF_A) {
		res = kmer.substr(0, posOfError + 1) + "A" + kmer.substr(posOfError + 1, kmer.size());
	} else if (error == ErrorType::DEL_OF_C) {
		res = kmer.substr(0, posOfError + 1) + "C" + kmer.substr(posOfError + 1, kmer.size());
	} else if (error == ErrorType::DEL_OF_G) {
		res = kmer.substr(0, posOfError + 1) + "G" + kmer.substr(posOfError + 1, kmer.size());
	} else if (error == ErrorType::DEL_OF_T) {
		res = kmer.substr(0, posOfError + 1) + "T" + kmer.substr(posOfError + 1, kmer.size());
	} else {
		throw std::runtime_error("Wrong error type given: " + errorTypeToString(error));
	}

	return res;
}

// TODO: This currently only accounts for single-base errors.
double KmerCounter::countKmerNoRCApproximate(const std::string &kmer,
		const std::shared_ptr<ErrorProfileUnit> &errorProfile) {
	/*if ((kmer.size() < 17) && (bufferApprox.find(kmer) != bufferApprox.end())) {
		return bufferApprox[kmer];
	}*/
	double countOriginal = sdsl::count(fm_index, kmer.begin(), kmer.end());
	double countTotal = 0;
	auto profileInformation = errorProfile->getKmerErrorProbabilities(kmer);
	double probCorrect = 0;
	for (size_t i = 0; i < profileInformation.size(); ++i) {
		for (auto kv : profileInformation[i]) {
			if (kv.first != ErrorType::CORRECT && kv.first != ErrorType::NODEL && kv.first != ErrorType::MULTIDEL) {
				std::string kmerAfterCorrection = kmerAfterError(kmer, kv.first, i);
				countTotal += exp(kv.second)
						* sdsl::count(fm_index, kmerAfterCorrection.begin(), kmerAfterCorrection.end());
			}
		}
		probCorrect += profileInformation[i][ErrorType::CORRECT] + profileInformation[i][ErrorType::NODEL];
	}
	countTotal += exp(probCorrect) * countOriginal;

	/*if (kmer.size() < 17 && countTotal >= 50) {
		bufferApprox[kmer] = countTotal;
	}*/
	return countTotal;
}

std::string KmerCounter::reverseComplementString(const std::string &sequence) {
	std::string rc = "";
	for (int i = sequence.size() - 1; i >= 0; i--) {
		assert(i >= 0);
		if (sequence[i] == 'A') {
			rc += "T";
		} else if (sequence[i] == 'T') {
			rc += "A";
		} else if (sequence[i] == 'C') {
			rc += "G";
		} else if (sequence[i] == 'G') {
			rc += "C";
		} else if (sequence[i] == 'N') {
			rc += "N";
		} else {
			throw std::runtime_error("Error:" + std::to_string(sequence[i]) + " at pos " + std::to_string(i) + " of " + std::to_string(sequence.size())+ " invalid.");
		}
	}
	return rc;
}
