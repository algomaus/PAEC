/*
 * ErrorCorrectionUnit.cpp
 *
 *  Created on: Feb 6, 2017
 *      Author: sarah
 */

#include "ErrorCorrectionUnit.h"

#include <functional>
#include <iostream>

#include <cereal/archives/binary.hpp>

#include "../ProducerConsumerPattern.hpp"
#include "../CorrectedRead.h"

ErrorCorrectionUnit::ErrorCorrectionUnit() {
}

ErrorCorrectionUnit::ErrorCorrectionUnit(ErrorCorrectionType type, ErrorProfileUnit &epu, KmerClassificationUnit &kcu, bool correctIndels) {
	if (type == ErrorCorrectionType::KMER_BASED) {
		correctRead = std::bind(correctRead_KmerBased, _1, std::ref(epu), std::ref(kcu), correctIndels);
	} else if (type == ErrorCorrectionType::NAIVE) {
		correctRead = std::bind(correctRead_Naive, _1, std::ref(epu), std::ref(kcu), correctIndels);
	} else {
		throw std::runtime_error("Unclear error correction type");
	}
}

void ErrorCorrectionUnit::addReadsFile(const std::string &filepath) {
	readFiles.push_back(filepath);
	outFilesCorrectedReads.push_back(std::ofstream(filepath + ".correctedReads.fastq"));
	outFilesCorrections.push_back(std::ofstream(filepath + ".corrections.txt"));
	iterators.push_back(std::make_unique<FASTQModifiedIterator>(filepath));

	cereal::BinaryOutputArchive oarchive(outFilesCorrections[outFilesCorrections.size() - 1]);
	oarchive(iterators[iterators.size() - 1]->numReadsLeft());
}

void ErrorCorrectionUnit::addReadsFile(const std::string &filepath, const std::string &outputPath) {
	readFiles.push_back(filepath);
	outFilesCorrectedReads.push_back(std::ofstream(outputPath + "correctedReads.fastq"));
	outFilesCorrections.push_back(std::ofstream(outputPath + "corrections.txt"));
	iterators.push_back(std::make_unique<FASTQModifiedIterator>(filepath));

	cereal::BinaryOutputArchive oarchive(outFilesCorrections[outFilesCorrections.size() - 1]);
	oarchive(iterators[iterators.size() - 1]->numReadsLeft());
}

void ErrorCorrectionUnit::correctReads() {
	for (size_t i = 0; i < readFiles.size(); ++i) {
		double minProgress = 0;
		while (iterators[i]->hasReadsLeft()) {
			FASTQRead fastqRead = iterators[i]->next();
			CorrectedRead cr = correctRead(fastqRead);
			if (cr.correctedRead.sequence.empty()) {
				throw std::runtime_error("The corrected read is empty!");
			}
			outFilesCorrectedReads[i] << cr.correctedRead << "\n";
			//if (correctedRead.corrections.size() > 0) {
			cereal::BinaryOutputArchive oarchive(outFilesCorrections[i]);
			oarchive(cr);
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

void ErrorCorrectionUnit::correctReadsMultithreaded() {
	outMtx = std::vector<std::mutex>(readFiles.size());
	auto fpProduce = std::bind(&ErrorCorrectionUnit::produceData, this, _1, _2);
	auto fpConsume = std::bind(&ErrorCorrectionUnit::consumeData, this, _1, _2);

	ProducerConsumerPattern<FASTQRead> pct(50, fpProduce, fpConsume);
	pct.run(readFiles.size(), consumersPerFile * readFiles.size());

	for (size_t i = 0; i < readFiles.size(); ++i) {
		outFilesCorrectedReads[i].close();
		outFilesCorrections[i].close();
	}
}

double ErrorCorrectionUnit::produceData(std::vector<FASTQRead> &buffer, size_t producerId) {
	if (iterators[producerId]->hasReadsLeft()) {
		std::vector<FASTQRead> fastqReads = iterators[producerId]->next(maxBufferSize);
		for (FASTQRead f : fastqReads) {
			buffer.push_back(f);
		}
	}

	return iterators[producerId]->progress();
}

void ErrorCorrectionUnit::consumeData(std::vector<FASTQRead> &buffer, size_t consumerId) {
	for (FASTQRead fastqRead : buffer) {
		CorrectedRead cr = correctRead(fastqRead);

		if (cr.correctedRead.sequence.empty()) {
			throw std::runtime_error("The corrected read is empty!");
		}

		std::string correctedReadString;

		std::stringstream ss;
		ss << cr.correctedRead;
		correctedReadString = ss.str();

		std::lock_guard<std::mutex> lck(outMtx[consumerId / consumersPerFile]);

		outFilesCorrectedReads[consumerId / consumersPerFile] << correctedReadString + "\n";
		//if (cra.corrections.size() > 0) {
		cereal::BinaryOutputArchive oarchive(outFilesCorrections[consumerId / consumersPerFile]);
		oarchive(cr);
		//}
	}
}
