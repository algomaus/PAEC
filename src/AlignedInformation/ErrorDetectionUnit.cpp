/*
 * ErrorDetectionUnit.cpp
 *
 *  Created on: Jan 31, 2017
 *      Author: sarah
 */

#include "ErrorDetectionUnit.h"

#include <cereal/archives/binary.hpp>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <cassert>
#include <functional>
#include <iostream>
#include <fstream>

#include "../FASTQRead.h"
#include "../ProducerConsumerPattern.hpp"
#include "CorrectedReadAligned.h"

using namespace std::placeholders;

ErrorDetectionUnit::ErrorDetectionUnit() {
}

void ErrorDetectionUnit::addAlignmentsFile(const std::string &alignmentFilePath) {
	alignmentsFiles.push_back(alignmentFilePath);
	outFilesCorrectedReads.push_back(std::ofstream(alignmentFilePath + ".trueReads.fasta"));
	outFilesCorrections.push_back(std::ofstream(alignmentFilePath + ".trueCorrections.txt", std::ios::binary));
	//BAMIterator it(alignmentFilePath, genome);
	iterators.push_back(std::make_unique<BAMIterator>(alignmentFilePath));

	cereal::BinaryOutputArchive oarchive(outFilesCorrections[outFilesCorrections.size() - 1]);
	oarchive(iterators[iterators.size() - 1]->getNumReadsTotalMapped());
}

void ErrorDetectionUnit::addAlignmentsFile(const std::string &alignmentFilePath, const std::string &outputPath) {
	alignmentsFiles.push_back(alignmentFilePath);
	outFilesCorrectedReads.push_back(std::ofstream(outputPath + "trueReads.fasta"));
	outFilesCorrections.push_back(std::ofstream(outputPath + "trueCorrections.txt", std::ios::binary));
	iterators.push_back(std::make_unique<BAMIterator>(alignmentFilePath));

	cereal::BinaryOutputArchive oarchive(outFilesCorrections[outFilesCorrections.size() - 1]);
	oarchive(iterators[iterators.size() - 1]->getNumReadsTotalMapped());
}

void ErrorDetectionUnit::correctReads(const seqan::Dna5String &reference) {
	genomePtr = std::make_shared<seqan::Dna5String>(reference);
	for (size_t i = 0; i < alignmentsFiles.size(); ++i) {
		double minProgress = 1;
		while (iterators[i]->hasReadsLeft()) {
			ReadWithAlignments alignedRead = iterators[i]->next();
			CorrectedReadAligned cra = alignedRead.retrieveCorrectedRead(reference);
			outFilesCorrectedReads[i] << cra.correctedRead << "\n";
			//if (cra.corrections.size() > 0) {

			cereal::BinaryOutputArchive oarchive(outFilesCorrections[i]);
			oarchive(cra);
			//}

			double progress = iterators[i]->progress();
			if (progress >= minProgress) {
				std::cout << progress << " \%" << std::endl;
				minProgress += 1;
			}
		}
		outFilesCorrectedReads[i].close();
		outFilesCorrections[i].close();
	}
}

void ErrorDetectionUnit::correctReadsMultithreaded(const seqan::Dna5String &reference) {
	outMtx = std::vector<std::mutex>(alignmentsFiles.size());
	auto fpProduce = std::bind(&ErrorDetectionUnit::produceData, this, _1, _2);
	auto fpConsume = std::bind(&ErrorDetectionUnit::consumeData, this, _1, _2);

	genomePtr = std::make_shared<seqan::Dna5String>(reference);
	ProducerConsumerPattern<ReadWithAlignments> pct(50, fpProduce, fpConsume);
	pct.run(alignmentsFiles.size(), consumersPerFile * alignmentsFiles.size());

	for (size_t i = 0; i < alignmentsFiles.size(); ++i) {
		outFilesCorrectedReads[i].close();
		outFilesCorrections[i].close();
	}
}

double ErrorDetectionUnit::produceData(std::vector<ReadWithAlignments> &buffer, size_t producerId) {
	if (iterators[producerId]->hasReadsLeft()) {
		std::vector<ReadWithAlignments> alignedReads = iterators[producerId]->next(maxBufferSize);
		for (ReadWithAlignments rwa : alignedReads) {
			buffer.push_back(rwa);
		}
	}

	return iterators[producerId]->progress();
}

void ErrorDetectionUnit::consumeData(std::vector<ReadWithAlignments> &buffer, size_t consumerId) {
	assert(genomePtr != NULL);
	for (ReadWithAlignments alignedRead : buffer) {
		CorrectedReadAligned cra = alignedRead.retrieveCorrectedRead(*genomePtr.get());

		std::string correctedReadString;
		std::string craString;

		std::stringstream ss;
		ss << cra.correctedRead;
		correctedReadString = ss.str();

		std::lock_guard<std::mutex> lck(outMtx[consumerId / consumersPerFile]);

		outFilesCorrectedReads[consumerId / consumersPerFile] << correctedReadString << "\n";
		//if (cra.corrections.size() > 0) {
		cereal::BinaryOutputArchive oarchive(outFilesCorrections[consumerId / consumersPerFile]);
		oarchive(cra);
		//}
	}
}
