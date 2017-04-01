/*
 * FeatureExtractor.hpp
 *
 *  Created on: Feb 8, 2017
 *      Author: sarah
 */

#pragma once

#include <memory>
#include "../../Correction.h"
#include "../../ErrorType.h"
#include "../../FASTQRead.h"
#include "../../KmerClassification/KmerClassificationUnit.h"
#include "../motif_analysis/MotifErrorProfile.h"

class FeatureExtractor {
public:
	virtual ~FeatureExtractor() {
	}
	FeatureExtractor(KmerClassificationUnit &kmerClassificationUnit, MotifErrorProfile &motifProfile, bool useKmerZScoresOrNot = false, bool useErrorsOnly = true) :
			kmerClassifier(kmerClassificationUnit), mep(motifProfile) {
		minKmerSize = kmerClassificationUnit.getMinKmerSize();
		useKmerZScores = useKmerZScoresOrNot;
		errorsOnly = useErrorsOnly;
		//lastKmerSize = minKmerSize;
	}
	virtual std::vector<double> getFeatureVector(const FASTQRead &read, size_t posInRead) = 0;
	virtual std::vector<double> getFeatureVectorNoQual(const std::string &sequence, size_t posInSequence) = 0;
	std::string getTrainingString(const FASTQRead &read, size_t posInRead, ErrorType type) {
		std::string data;
		std::vector<double> features = getFeatureVector(read, posInRead);
		for (size_t i = 0; i < features.size(); ++i) {
			data += std::to_string(features[i]) + ";";
		}
		data += std::to_string(errorTypeToNumber(type));
		return data;
	}
	std::string getTrainingStringNoQual(const std::string &sequence, size_t posInSequence, ErrorType type) {
		std::string data;
		std::vector<double> features = getFeatureVectorNoQual(sequence, posInSequence);
		for (size_t i = 0; i < features.size(); ++i) {
			data += std::to_string(features[i]) + ";";
		}
		data += std::to_string(errorTypeToNumber(type));
		return data;
	}
	virtual std::string getFeatureNames() {
		return featureNames;
	}
	virtual std::string getFeatureNamesNoQual() {
		return featureNamesNoQual;
	}
	virtual std::vector<std::string> getFeatureNamesVector() {
		return featureNamesVector;
	}
	virtual std::vector<std::string> getFeatureNamesNoQualVector() {
		return featureNamesNoQualVector;
	}
	virtual std::vector<ErrorType> getClasses() {
		return classes;
	}
protected:
	// TODO: Make this way faster! Maybe by using some kind of binary search?
	double kmerZScoreExtract(const std::string &sequence, size_t posInRead, const std::string &middleAs,
			bool leftPossibleOrig, bool rightPossibleOrig) {
		bool leftPossible = leftPossibleOrig;
		bool rightPossible = rightPossibleOrig;

		double zScore = 0;
		std::string kmer = middleAs;
		size_t offsetLeft = 1;
		size_t offsetRight = 1;
		bool goLeft = leftPossible;
		size_t maxKmerSize = sequence.size();

		/*
		// trying to speed up things: start
		if (lastKmerSize > minKmerSize) {
			while ((leftPossible || rightPossible) && kmer.size() < lastKmerSize - 2) {
				if (goLeft) {
					if (offsetLeft > posInRead) {
						leftPossible = false;
						if (rightPossible) {
							goLeft = false;
						}
					} else {
						kmer = sequence[posInRead - offsetLeft] + kmer;
						offsetLeft++;
					}
				} else {
					if (posInRead + offsetRight >= sequence.size()) {
						rightPossible = false;
						if (leftPossible) {
							goLeft = true;
						}
					} else {
						kmer = sequence[posInRead + offsetRight] + kmer;
						offsetRight++;
					}
				}
			}
			std::cout << "last k-mer size: " << lastKmerSize << "\n";
			if (kmerClassifier.classifyKmer(kmer) != KmerType::REPEAT) { // speed up did not work -> restart
				leftPossible = leftPossibleOrig;
				rightPossible = rightPossibleOrig;
				zScore = 0;
				kmer = middleAs;
				offsetLeft = 1;
				offsetRight = 1;
				goLeft = leftPossible;
				std::cout << "Speedup failed :-(";
			}
		}
		// trying to speed up things: end
		*/

		while ((leftPossible || rightPossible) && kmer.size() < maxKmerSize) {
			if (goLeft) {
				if (offsetLeft > posInRead) {
					leftPossible = false;
					if (rightPossible) {
						goLeft = false;
					}
				} else {
					kmer = sequence[posInRead - offsetLeft] + kmer;
					if (kmer.size() >= minKmerSize && kmer.size() % 2 == 1) {
						if (kmerClassifier.classifyKmer(kmer) != KmerType::REPEAT) {
							zScore = kmerClassifier.kmerZScore(kmer);
							//lastKmerSize = kmer.size();
							break;
						}
					}
					offsetLeft++;
				}
			} else {
				if (posInRead + offsetRight >= sequence.size()) {
					rightPossible = false;
					if (leftPossible) {
						goLeft = true;
					}
				} else {
					kmer = sequence[posInRead + offsetRight] + kmer;
					if (kmer.size() >= minKmerSize && kmer.size() % 2 == 1) {
						if (kmerClassifier.classifyKmer(kmer) != KmerType::REPEAT) {
							zScore = kmerClassifier.kmerZScore(kmer);
							//lastKmerSize = kmer.size();
							break;
						}
					}
					offsetRight++;
				}
			}
		}
		return zScore;
	}

	double motifZScoreExtract(const std::string &sequence, size_t posInRead, ErrorType type) {
		return mep.findMostSignificantZScore(type, sequence, posInRead);
	}

	std::string featureNames, featureNamesNoQual;
	std::vector<std::string> featureNamesVector, featureNamesNoQualVector;
	std::vector<ErrorType> classes;
	KmerClassificationUnit &kmerClassifier;
	MotifErrorProfile &mep;
	size_t minKmerSize;

	bool useKmerZScores;
	bool errorsOnly;

	//size_t lastKmerSize;
};

