/*
 * CoverageBiasUnit.h
 *
 *  Created on: Jan 31, 2017
 *      Author: sarah
 */

#pragma once

#include <stddef.h>
#include <string>
#include <cmath>
#include <memory>
#include <seqan/sequence.h>
#include <cereal/types/vector.hpp>

#include "../KmerClassification/KmerCounter.h"
#include "PUSM.h"
#include "../external/IntervalTree.h"


enum CoverageBiasType {
	IGNORE_BIAS = 0, MEDIAN_BIAS_ALIGNMENT = 1, MEDIAN_BIAS_REFERENCE = 2, MEDIAN_BIAS_READS_ONLY = 3
};


class CoverageBiasUnit {
public:
	CoverageBiasUnit();
	CoverageBiasUnit(CoverageBiasType coverageBiasType, size_t estimatedGenomeSize, PerfectUniformSequencingModel &pusm);
	double getBias(const std::string &kmer);
	void loadBias(const std::string &filepath);
	void storeBias(const std::string &filepath);
	void learnBiasFromReferenceAlignment(const seqan::Dna5String &referenceGenome, KmerCounter &referenceCounter,
			const std::string &alignmentsFilename, const std::unordered_map<size_t, size_t> &readLengths);
	void learnBiasFromReferenceMatches(const seqan::Dna5String &referenceGenome, KmerCounter &referenceCounter,
			KmerCounter &readsCounter, const std::unordered_map<size_t, size_t> &readLengths);
	void learnBiasFromReadsOnly(const std::string &readsFilePath, KmerCounter &readsCounter,
			const std::unordered_map<size_t, size_t> &readLengths);
	void printBias();
	void plotBias(const std::string &filename);

	size_t getMinKmerSize() {
		return minKmerSize;
	}

	std::vector<double> getBiases() {
		return biases;
	}

	double computeGCContent(const std::string &sequence);

	template<class Archive>
		void serialize(Archive & archive) {
			archive(minKmerSize, gc_step, biases); // serialize things by passing them to the archive
	}
protected:
	void fixEmptyBiases();

	size_t minKmerSize;
	size_t genomeSize;
	std::vector<double> biases;
	double gc_step;

	CoverageBiasType covBiasType;

	std::shared_ptr<PerfectUniformSequencingModel> pusm;
};
