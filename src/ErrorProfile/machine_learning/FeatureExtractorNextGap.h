/*
 * FeatureExtractorNextGap.h
 *
 *  Created on: Feb 8, 2017
 *      Author: sarah
 */

#pragma once

#include "FeatureExtractor.hpp"

class FeatureExtractorNextGap : public FeatureExtractor {
public:
	FeatureExtractorNextGap(KmerClassificationUnit &kmerClassifier, MotifErrorProfile &motifProfile, bool useKmerZScores = false, bool errorsOnly = true);
	virtual std::vector<double> getFeatureVector(const FASTQRead &read, size_t posInRead);
	virtual std::vector<double> getFeatureVectorNoQual(const std::string &sequence, size_t posInSequence);
};
