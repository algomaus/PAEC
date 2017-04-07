/*
 * ErrorDetectionUnit.h
 *
 *  Created on: Jan 31, 2017
 *      Author: sarah
 */

#pragma once

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <stddef.h>
#include <fstream>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include "../ErrorCorrectionEvaluation.h"
#include "ReadWithAlignments.h"

#include "BAMIterator.h"

#include "../ErrorProfile/ErrorProfileUnit.hpp"

/*
 * "Corrects" the reads by using the mappings to the known reference genome.
 */

class ErrorDetectionUnit {
public:
	ErrorDetectionUnit();
	ErrorDetectionUnit(ErrorCorrectionEvaluation &ece);

	void addAlignmentsFile(const std::string &alignmentFilePath);
	void addAlignmentsFile(const std::string &alignmentFilePath, const std::string &outputPath);
	void correctReads(const seqan::Dna5String &reference);
	void correctReadsMultithreaded(const seqan::Dna5String &reference);

	void addObserver(ErrorProfileUnit& epuObs);
private:
	double produceData(std::vector<ReadWithAlignments> &buffer, size_t producerId);
	void consumeData(std::vector<ReadWithAlignments> &buffer, size_t consumerId);

	std::vector<std::string> alignmentsFiles;
	std::vector<std::ofstream> outFilesCorrectedReads;
	//std::vector<std::ofstream> outFilesCorrections;
	std::vector<std::mutex> outMtx;
	std::vector<std::unique_ptr<BAMIterator> > iterators;
	std::shared_ptr<seqan::Dna5String> genomePtr;

	std::vector<ErrorProfileUnit*> observers;
	ErrorCorrectionEvaluation* ecEval;

	size_t consumersPerFile = 3;
	size_t maxBufferSize = consumersPerFile * 2;
};
