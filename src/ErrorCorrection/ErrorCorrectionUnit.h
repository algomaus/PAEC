/*
 * ErrorCorrectionUnit.h
 *
 *  Created on: Jan 31, 2017
 *      Author: sarah
 */

#pragma once

#include <stddef.h>
#include <fstream>
#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include "../CorrectedRead.h"
#include "../ErrorCorrectionEvaluation.h"
#include "../FASTQRead.h"
#include "../ErrorProfile/ErrorProfileUnit.hpp"
#include "../FASTQModifiedIterator.h"
#include "../KmerClassification/KmerClassificationUnit.h"
#include "ErrorCorrectionAlgorithms.h"

using namespace std::placeholders;

enum ErrorCorrectionType {
	NAIVE = 0, KMER_BASED = 1, KMER_IMPROVED = 2
};

/*
 * TODO: This class has lots of code duplication with ErrorDetectionUnit.h
 */

class ErrorCorrectionUnit {
public:
	ErrorCorrectionUnit();
	ErrorCorrectionUnit(ErrorCorrectionType type, ErrorProfileUnit &epu, KmerClassificationUnit &kcu, bool correctIndels, ErrorCorrectionEvaluation &ece);
	void addReadsFile(const std::string &filepath);
	void addReadsFile(const std::string &filepath, const std::string &outputPath);
	void correctReads(); // write the results to filepath + "_precorrected.txt"
	void correctReadsMultithreaded();

	void addObserver(ErrorProfileUnit& epuObs);
private:
	double produceData(std::vector<FASTQRead> &buffer, size_t producerId);
	void consumeData(std::vector<FASTQRead> &buffer, size_t consumerId);

	std::vector<std::string> readFiles;
	std::vector<std::ofstream> outFilesCorrectedReads;
	//std::vector<std::ofstream> outFilesCorrections;
	std::vector<std::mutex> outMtx;
	std::vector<std::unique_ptr<FASTQModifiedIterator> > iterators;
	size_t consumersPerFile = 3;
	size_t maxBufferSize = consumersPerFile * 20;

	std::vector<ErrorProfileUnit*> observers;

	ErrorCorrectionEvaluation* ecEval;

	std::function<CorrectedRead(FASTQRead)> correctRead;
};
