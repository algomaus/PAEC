/*
 * FeatureExtractorNextGap.cpp
 *
 *  Created on: Feb 8, 2017
 *      Author: sarah
 */

#include "FeatureExtractorNextGap.h"

FeatureExtractorNextGap::FeatureExtractorNextGap(KmerClassificationUnit &kmerClassifier,
		MotifErrorProfile &motifProfile, bool useKmerZScoresOrNot, bool useOnlyErrors) :
		FeatureExtractor(kmerClassifier, motifProfile, useKmerZScoresOrNot, useOnlyErrors) {
	for (ErrorType type : errorTypesNextGap()) {
		if (useOnlyErrors && type == ErrorType::NODEL) continue;
		classes.push_back(type);
	}

	if (useKmerZScores) {
		featureNames =
				"current base;position in read;read length;quality score;k-merZ left '*_';k-merZ middle '*_';k-merZ right '*_';k-merZ left '*A';k-merZ middle '*A';k-merZ right '*A';k-merZ left '*C';k-merZ middle '*C';k-merZ right '*C';k-merZ left '*G';k-merZ middle '*G';k-merZ right '*G';k-merZ left '*T';k-merZ middle '*T';k-merZ right '*T';MULTIDEL motifZ;DEL_OF_A motifZ;DEL_OF_C motifZ;DEL_OF_G motifZ;DEL_OF_T motifZ;type";
		featureNamesNoQual =
				"current base;position in read;read length;k-merZ left '*_';k-merZ middle '*_';k-merZ right '*_';k-merZ left '*A';k-merZ middle '*A';k-merZ right '*A';k-merZ left '*C';k-merZ middle '*C';k-merZ right '*C';k-merZ left '*G';k-merZ middle '*G';k-merZ right '*G';k-merZ left '*T';k-merZ middle '*T';k-merZ right '*T';MULTIDEL motifZ;DEL_OF_A motifZ;DEL_OF_C motifZ;DEL_OF_G motifZ;DEL_OF_T motifZ;type";
	} else {
		featureNames =
				"current base;position in read;read length;quality score;MULTIDEL motifZ;DEL_OF_A motifZ;DEL_OF_C motifZ;DEL_OF_G motifZ;DEL_OF_T motifZ;type";
		featureNamesNoQual =
				"current base;position in read;read length;MULTIDEL motifZ;DEL_OF_A motifZ;DEL_OF_C motifZ;DEL_OF_G motifZ;DEL_OF_T motifZ;type";
	}

	featureNamesVector.push_back("current base");
	featureNamesVector.push_back("position in read");
	featureNamesVector.push_back("read length");
	featureNamesVector.push_back("quality score");
	if (useKmerZScores) {
		featureNamesVector.push_back("k-merZ left '*_'");
		featureNamesVector.push_back("k-merZ middle '*_'");
		featureNamesVector.push_back("k-merZ right '*_'");
		featureNamesVector.push_back("k-merZ left '*A'");
		featureNamesVector.push_back("k-merZ middle '*A'");
		featureNamesVector.push_back("k-merZ right '*A'");
		featureNamesVector.push_back("k-merZ left '*C'");
		featureNamesVector.push_back("k-merZ middle '*C'");
		featureNamesVector.push_back("k-merZ right '*C'");
		featureNamesVector.push_back("k-merZ left '*G'");
		featureNamesVector.push_back("k-merZ middle '*G'");
		featureNamesVector.push_back("k-merZ right '*G'");
		featureNamesVector.push_back("k-merZ left '*T'");
		featureNamesVector.push_back("k-merZ middle '*T'");
		featureNamesVector.push_back("k-merZ right '*T'");
	}
	featureNamesVector.push_back("MULTIDEL motifZ");
	featureNamesVector.push_back("DEL_OF_A motifZ");
	featureNamesVector.push_back("DEL_OF_C motifZ");
	featureNamesVector.push_back("DEL_OF_G motifZ");
	featureNamesVector.push_back("DEL_OF_T motifZ");

	featureNamesNoQualVector.push_back("current base");
	featureNamesNoQualVector.push_back("position in read");
	featureNamesNoQualVector.push_back("read length");
	if (useKmerZScores) {
		featureNamesNoQualVector.push_back("k-merZ left '*_'");
		featureNamesNoQualVector.push_back("k-merZ middle '*_'");
		featureNamesNoQualVector.push_back("k-merZ right '*_'");
		featureNamesNoQualVector.push_back("k-merZ left '*A'");
		featureNamesNoQualVector.push_back("k-merZ middle '*A'");
		featureNamesNoQualVector.push_back("k-merZ right '*A'");
		featureNamesNoQualVector.push_back("k-merZ left '*C'");
		featureNamesNoQualVector.push_back("k-merZ middle '*C'");
		featureNamesNoQualVector.push_back("k-merZ right '*C'");
		featureNamesNoQualVector.push_back("k-merZ left '*G'");
		featureNamesNoQualVector.push_back("k-merZ middle '*G'");
		featureNamesNoQualVector.push_back("k-merZ right '*G'");
		featureNamesNoQualVector.push_back("k-merZ left '*T'");
		featureNamesNoQualVector.push_back("k-merZ middle '*T'");
		featureNamesNoQualVector.push_back("k-merZ right '*T'");
	}
	featureNamesNoQualVector.push_back("MULTIDEL motifZ");
	featureNamesNoQualVector.push_back("DEL_OF_A motifZ");
	featureNamesNoQualVector.push_back("DEL_OF_C motifZ");
	featureNamesNoQualVector.push_back("DEL_OF_G motifZ");
	featureNamesNoQualVector.push_back("DEL_OF_T motifZ");
}

std::vector<double> FeatureExtractorNextGap::getFeatureVector(const FASTQRead &read, size_t posInRead) {
	std::vector<double> vec;
	vec.push_back((double) read.sequence[posInRead]);
	vec.push_back((double) posInRead);
	vec.push_back((double) read.sequence.size());
	vec.push_back((double) read.quality[posInRead]);

	std::string actBase;
	actBase += read.sequence[posInRead];

	double mZ1, mZ2, mZ3, mZ4, mZ5;

	if (useKmerZScores) {
		std::vector<double> kmerZScores(15);

		std::vector<std::string> middle = { actBase + "", actBase + "A", actBase + "C", actBase + "G", actBase + "T" };

		for (size_t i = 0; i < middle.size(); ++i) {
			kmerZScores[3 * i] = kmerZScoreExtract(read.sequence, posInRead, middle[i], true, false);
			kmerZScores[3 * i + 1] = kmerZScoreExtract(read.sequence, posInRead, middle[i], true, true);
			kmerZScores[3 * i + 2] = kmerZScoreExtract(read.sequence, posInRead, middle[i], false, true);
		}

		for (size_t i = 0; i < kmerZScores.size(); ++i) {
			vec.push_back(kmerZScores[i]);
		}
	}

	mZ1 = motifZScoreExtract(read.sequence, posInRead, ErrorType::MULTIDEL);
	mZ2 = motifZScoreExtract(read.sequence, posInRead, ErrorType::DEL_OF_A);
	mZ3 = motifZScoreExtract(read.sequence, posInRead, ErrorType::DEL_OF_C);
	mZ4 = motifZScoreExtract(read.sequence, posInRead, ErrorType::DEL_OF_G);
	mZ5 = motifZScoreExtract(read.sequence, posInRead, ErrorType::DEL_OF_T);

	vec.push_back(mZ1);
	vec.push_back(mZ2);
	vec.push_back(mZ3);
	vec.push_back(mZ4);
	vec.push_back(mZ5);

	return vec;
}

std::vector<double> FeatureExtractorNextGap::getFeatureVectorNoQual(const std::string &sequence, size_t posInSequence) {
	std::vector<double> vec;
	vec.push_back((double) sequence[posInSequence]);
	vec.push_back((double) posInSequence);
	vec.push_back((double) sequence.size());

	std::string actBase;
	actBase += sequence[posInSequence];

	if (useKmerZScores) {
		std::vector<std::string> middle = { actBase + "", actBase + "A", actBase + "C", actBase + "G", actBase + "T" };
		for (size_t i = 0; i < middle.size(); ++i) {
			vec.push_back(kmerZScoreExtract(sequence, posInSequence, middle[i], true, false));
			vec.push_back(kmerZScoreExtract(sequence, posInSequence, middle[i], true, true));
			vec.push_back(kmerZScoreExtract(sequence, posInSequence, middle[i], false, true));
		}
	}

	vec.push_back(motifZScoreExtract(sequence, posInSequence, ErrorType::MULTIDEL));
	vec.push_back(motifZScoreExtract(sequence, posInSequence, ErrorType::DEL_OF_A));
	vec.push_back(motifZScoreExtract(sequence, posInSequence, ErrorType::DEL_OF_C));
	vec.push_back(motifZScoreExtract(sequence, posInSequence, ErrorType::DEL_OF_G));
	vec.push_back(motifZScoreExtract(sequence, posInSequence, ErrorType::DEL_OF_T));
	return vec;
}
