/*
 * KmerClassificationUnit.h
 *
 *  Created on: Jan 31, 2017
 *      Author: sarah
 */

#pragma once

#include <seqan/sequence.h>
#include <stddef.h>
#include <iostream>
#include <memory>
#include <string>
#include <python2.7/Python.h>

#include "../AlignedInformation/Dataset.hpp"
#include "KmerType.h"

#include "../CoverageBias/CoverageBiasUnit.h"
#include "KmerCounter.h"
#include "../CoverageBias/PUSM.h"

enum KmerClassificationType {
	CLASSIFICATION_STATISTICAL = 0, CLASSIFICATION_NAIVE = 1, CLASSIFICATION_MACHINE_LEARNING = 2
};

class KmerClassificationUnit {
public:
	KmerClassificationUnit(KmerCounter &kmerCounter, CoverageBiasUnit &biasUnitRef,
			PerfectUniformSequencingModel &pusmRef, KmerClassificationType type);
	~KmerClassificationUnit();
	void trainClassifier(Dataset &ds, KmerCounter &genomeCounter);
	KmerType classifyKmer(const std::string &kmer);
	KmerType classifyZScore(double zScore);
	double kmerZScore(const std::string &kmer);
	size_t getMinKmerSize();
	void storeClassifier(const std::string &filename);
	void loadClassifier(const std::string &filename);
	void clearCache();
private:
	void extractTrainingData(Dataset &ds, size_t k, std::ofstream &outfile);
	void extractTrainingDataFromReference(const seqan::Dna5String &referenceGenome, size_t k,
			KmerCounter &referenceCounter, std::ofstream &outfile);
	void writeTrainingString(double zScore, double gc, size_t k, double countObserved, double countBiasCorrected,
			double countExpectedPusm, KmerType type, std::ofstream &outfile);
	KmerCounter &counter;
	CoverageBiasUnit &biasUnit;
	PerfectUniformSequencingModel &pusm;
	KmerClassificationType classificationType;
	std::string featureNames;
	std::string outputPath;

	std::unordered_map<std::string, KmerType> cachedClassifications;

	PyObject* mlClassifier;
};
