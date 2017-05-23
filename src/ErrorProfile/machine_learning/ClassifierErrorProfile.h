/*
 * ClassifierErrorProfile.h
 *
 *  Created on: Jan 31, 2017
 *      Author: sarah
 */

#pragma once

#include <stddef.h>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include "../../external/cereal/types/unordered_map.hpp"
#include "../../external/cereal/types/utility.hpp"
#include <python2.7/Python.h>

#include "../../ErrorType.h"
#include "../../UtitilyFunctions.hpp"
#include "../ErrorProfileUnit.hpp"
#include "FeatureExtractorCurrentBase.h"
#include "FeatureExtractorNextGap.h"
#include "../../KmerClassification/KmerClassificationUnit.h"
#include "../../CorrectedRead.h"
#include "../motif_analysis/MotifErrorProfile.h"

class ClassifierErrorProfile: public ErrorProfileUnit {
public:
	ClassifierErrorProfile(const std::string &plotPath, CoverageBiasType coverageBiasType,
			KmerClassificationUnit &kmerClassifier, MotifErrorProfile &motifProfile, bool useQualityScores, bool useKmerZScores = false, bool useErrorsOnly = true);
	~ClassifierErrorProfile();
	virtual std::unordered_map<ErrorType, double> getErrorProbabilities(const FASTQRead &read, size_t positionInRead);
	virtual std::unordered_map<ErrorType, double> getKmerErrorProbabilities(const std::string &kmer,
			size_t positionInKmer);
	virtual void loadErrorProfile(const std::string &filepath, KmerCounter &counter);
	virtual void storeErrorProfile(const std::string &filepath);
	virtual void plotErrorProfile();

	virtual void reset();
	virtual void check(const CorrectedRead &corrRead, double acceptProb = 1.0);
	virtual void checkAligned(const CorrectedReadAligned &corrRead, double acceptProb = 1.0);

	void setOverallErrorRates(double errorRateBase, double errorRateGap);

	virtual void finalize();

	template<class Archive>
	void serialize(Archive & archive) {
		archive(clsfyCurrentBase, clsfyNextGap, trainCurrentBase, trainNextGap, overallErrorRateCurrentBase,
				overallErrorRateNextGap, finalized, useKmerZScores); // serialize things by passing them to the archive
	}
protected:
	ClassifierErrorProfile();
	virtual std::unordered_map<ErrorType, double> getErrorProbabilitiesFinalized(const FASTQRead &read,
			size_t positionInRead);
	virtual std::unordered_map<ErrorType, double> getErrorProbabilitiesFinalized(const std::string &kmer,
			size_t positionInKmer);
private:
	void processCorrection(const FASTQRead &read, size_t posInRead, ErrorType type, std::vector<bool> &nodel,
			std::vector<bool> &correct);
	void processCorrect(const FASTQRead &read, size_t posInRead);
	void processNodel(const FASTQRead &read, size_t posInRead);

	PyObject* classifierCurrentBase;
	PyObject* classifierNextGap;
	std::unique_ptr<FeatureExtractorCurrentBase> feCurrentBase;
	std::unique_ptr<FeatureExtractorNextGap> feNextGap;

	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution;

	std::string clsfyCurrentBase; // path to the trained sklearn classifiers
	std::string clsfyNextGap; // path to the trained sklearn classifiers
	std::string trainCurrentBase; // path to the csv file containing the training data
	std::string trainNextGap; // path to the csv file containing the training data
	double overallErrorRateCurrentBase;
	double overallErrorRateNextGap;
	bool finalized;
	bool useQual;
	bool useKmerZScores;
	bool errorsOnly;

	bool alreadyHasTrainingData;
};
