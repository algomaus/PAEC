/*
 * FeatureExtractorCurrentBase.cpp
 *
 *  Created on: Feb 8, 2017
 *      Author: sarah
 */

#include "FeatureExtractorCurrentBase.h"

FeatureExtractorCurrentBase::FeatureExtractorCurrentBase(KmerClassificationUnit &kmerClassifier,
		MotifErrorProfile &motifProfile, bool useKmerZScoresOrNot, bool useOnlyErrors) :
		FeatureExtractor(kmerClassifier, motifProfile, useKmerZScoresOrNot, useOnlyErrors) {
	for (ErrorType type : errorTypesCurrentBase()) {
		if (useOnlyErrors && type == ErrorType::CORRECT) continue;
		classes.push_back(type);
	}

	if (useKmerZScores) {
		featureNames =
				"current base;position in read;read length;quality score;k-merZ left '_';k-merZ middle '_';k-merZ right '_';k-merZ left 'A';k-merZ middle 'A';k-merZ right 'A';k-merZ left 'C';k-merZ middle 'C';k-merZ right 'C';k-merZ left 'G';k-merZ middle 'G';k-merZ right 'G';k-merZ left 'T';k-merZ middle 'T';k-merZ right 'T';INSERTION motifZ;SUB_FROM_A motifZ;SUB_FROM_C motifZ;SUB_FROM_G motifZ;SUB_FROM_T motifZ;type";
		featureNamesNoQual =
				"current base;position in read;read length;k-merZ left '_';k-merZ middle '_';k-merZ right '_';k-merZ left 'A';k-merZ middle 'A';k-merZ right 'A';k-merZ left 'C';k-merZ middle 'C';k-merZ right 'C';k-merZ left 'G';k-merZ middle 'G';k-merZ right 'G';k-merZ left 'T';k-merZ middle 'T';k-merZ right 'T';INSERTION motifZ;SUB_FROM_A motifZ;SUB_FROM_C motifZ;SUB_FROM_G motifZ;SUB_FROM_T motifZ;type";
	} else {
		featureNames =
				"current base;position in read;read length;quality score;INSERTION motifZ;SUB_FROM_A motifZ;SUB_FROM_C motifZ;SUB_FROM_G motifZ;SUB_FROM_T motifZ;type";
		featureNamesNoQual =
				"current base;position in read;read length;INSERTION motifZ;SUB_FROM_A motifZ;SUB_FROM_C motifZ;SUB_FROM_G motifZ;SUB_FROM_T motifZ;type";
	}

	featureNamesVector.push_back("current base");
	featureNamesVector.push_back("position in read");
	featureNamesVector.push_back("read length");
	featureNamesVector.push_back("quality score");
	if (useKmerZScores) {
		featureNamesVector.push_back("k-merZ left '_'");
		featureNamesVector.push_back("k-merZ middle '_'");
		featureNamesVector.push_back("k-merZ right '_'");
		featureNamesVector.push_back("k-merZ left 'A'");
		featureNamesVector.push_back("k-merZ middle 'A'");
		featureNamesVector.push_back("k-merZ right 'A'");
		featureNamesVector.push_back("k-merZ left 'C'");
		featureNamesVector.push_back("k-merZ middle 'C'");
		featureNamesVector.push_back("k-merZ right 'C'");
		featureNamesVector.push_back("k-merZ left 'G'");
		featureNamesVector.push_back("k-merZ middle 'G'");
		featureNamesVector.push_back("k-merZ right 'G'");
		featureNamesVector.push_back("k-merZ left 'T'");
		featureNamesVector.push_back("k-merZ middle 'T'");
		featureNamesVector.push_back("k-merZ right 'T'");
	}
	featureNamesVector.push_back("INSERTION motifZ");
	featureNamesVector.push_back("SUB_FROM_A motifZ");
	featureNamesVector.push_back("SUB_FROM_C motifZ");
	featureNamesVector.push_back("SUB_FROM_G motifZ");
	featureNamesVector.push_back("SUB_FROM_T motifZ");

	featureNamesNoQualVector.push_back("current base");
	featureNamesNoQualVector.push_back("position in read");
	featureNamesNoQualVector.push_back("read length");
	if (useKmerZScores) {
		featureNamesNoQualVector.push_back("k-merZ left '_'");
		featureNamesNoQualVector.push_back("k-merZ middle '_'");
		featureNamesNoQualVector.push_back("k-merZ right '_'");
		featureNamesNoQualVector.push_back("k-merZ left 'A'");
		featureNamesNoQualVector.push_back("k-merZ middle 'A'");
		featureNamesNoQualVector.push_back("k-merZ right 'A'");
		featureNamesNoQualVector.push_back("k-merZ left 'C'");
		featureNamesNoQualVector.push_back("k-merZ middle 'C'");
		featureNamesNoQualVector.push_back("k-merZ right 'C'");
		featureNamesNoQualVector.push_back("k-merZ left 'G'");
		featureNamesNoQualVector.push_back("k-merZ middle 'G'");
		featureNamesNoQualVector.push_back("k-merZ right 'G'");
		featureNamesNoQualVector.push_back("k-merZ left 'T'");
		featureNamesNoQualVector.push_back("k-merZ middle 'T'");
		featureNamesNoQualVector.push_back("k-merZ right 'T'");
	}
	featureNamesNoQualVector.push_back("INSERTION motifZ");
	featureNamesNoQualVector.push_back("SUB_FROM_A motifZ");
	featureNamesNoQualVector.push_back("SUB_FROM_C motifZ");
	featureNamesNoQualVector.push_back("SUB_FROM_G motifZ");
	featureNamesNoQualVector.push_back("SUB_FROM_T motifZ");
}

std::vector<double> FeatureExtractorCurrentBase::getFeatureVector(const FASTQRead &read, size_t posInRead) {
	std::vector<double> vec;
	vec.push_back((double) read.sequence[posInRead]);
	vec.push_back((double) posInRead);
	vec.push_back((double) read.sequence.size());
	vec.push_back((double) read.quality[posInRead]);

	std::vector<std::string> middle = { "", "A", "C", "G", "T" };

	double mZ1, mZ2, mZ3, mZ4, mZ5;

	if (useKmerZScores) {
		std::vector<double> kmerZScores(15);
		for (size_t i = 0; i < middle.size(); ++i) {
			kmerZScores[3 * i] = kmerZScoreExtract(read.sequence, posInRead, middle[i], true, false);
			kmerZScores[3 * i + 1] = kmerZScoreExtract(read.sequence, posInRead, middle[i], true, true);
			kmerZScores[3 * i + 2] = kmerZScoreExtract(read.sequence, posInRead, middle[i], false, true);
		}
		for (size_t i = 0; i < kmerZScores.size(); ++i) {
			vec.push_back(kmerZScores[i]);
		}
	}

	mZ1 = motifZScoreExtract(read.sequence, posInRead, ErrorType::INSERTION);
	mZ2 = motifZScoreExtract(read.sequence, posInRead, ErrorType::SUB_FROM_A);
	mZ3 = motifZScoreExtract(read.sequence, posInRead, ErrorType::SUB_FROM_C);
	mZ4 = motifZScoreExtract(read.sequence, posInRead, ErrorType::SUB_FROM_G);
	mZ5 = motifZScoreExtract(read.sequence, posInRead, ErrorType::SUB_FROM_T);

	vec.push_back(mZ1);
	vec.push_back(mZ2);
	vec.push_back(mZ3);
	vec.push_back(mZ4);
	vec.push_back(mZ5);
	return vec;
}

std::vector<double> FeatureExtractorCurrentBase::getFeatureVectorNoQual(const std::string &sequence,
		size_t posInSequence) {
	std::vector<double> vec;
	vec.push_back((double) sequence[posInSequence]);
	vec.push_back((double) posInSequence);
	vec.push_back((double) sequence.size());

	if (useKmerZScores) {
		std::vector<std::string> middle = { "", "A", "C", "G", "T" };
		for (size_t i = 0; i < middle.size(); ++i) {
			vec.push_back(kmerZScoreExtract(sequence, posInSequence, middle[i], true, false));
			vec.push_back(kmerZScoreExtract(sequence, posInSequence, middle[i], true, true));
			vec.push_back(kmerZScoreExtract(sequence, posInSequence, middle[i], false, true));
		}
	}
	vec.push_back(motifZScoreExtract(sequence, posInSequence, ErrorType::INSERTION));
	vec.push_back(motifZScoreExtract(sequence, posInSequence, ErrorType::SUB_FROM_A));
	vec.push_back(motifZScoreExtract(sequence, posInSequence, ErrorType::SUB_FROM_C));
	vec.push_back(motifZScoreExtract(sequence, posInSequence, ErrorType::SUB_FROM_G));
	vec.push_back(motifZScoreExtract(sequence, posInSequence, ErrorType::SUB_FROM_T));
	return vec;
}
