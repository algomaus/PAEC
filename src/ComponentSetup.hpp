/*
 * ComponentSetup.hpp
 *
 *  Created on: Mar 15, 2017
 *      Author: sarah
 */

#pragma once

#include <functional>
#include <memory>

#include "AlignedInformation/Dataset.hpp"
#include "AlignedInformation/ErrorDetectionUnit.h"
#include "CorrectedRead.h"
#include "CoverageBias/PUSM.h"
#include "ErrorCorrection/ErrorCorrectionUnit.h"
#include "ErrorProfile/machine_learning/ClassifierErrorProfile.h"
#include "ErrorProfile/motif_analysis/MotifErrorProfile.h"
#include "ErrorProfile/OverallErrorProfile.h"
#include "FASTQRead.h"
#include "KmerClassification/KmerClassificationUnit.h"
#include "KmerClassification/KmerCounter.h"

using namespace std::placeholders;

enum ErrorProfileType {
	OVERALL_STATS_ONLY = 0, MOTIF_STATS_ONLY = 1, MACHINE_LEARNING = 2
};

class ComponentSetup {
public:
	ComponentSetup(Dataset &ds, CoverageBiasType coverageBiasType, KmerClassificationType classificationType,
			ErrorProfileType errorProfileType, ErrorCorrectionType correctionType, bool correctIndels) :
			dataset(ds), counterReads(ds.readsOnlyFileName), counterReference(ds.referenceFileName), pusm(
					ds.readLengthsPtr, ds.genomeSize, ds.genomeType), covBias(coverageBiasType, ds.genomeSize, pusm), kmerClassifier(
					counterReads, covBias, pusm, classificationType), epuMotif(counterReads), epuClassify(ds.plotPath, coverageBiasType,
					kmerClassifier, epuMotif, ds.hasQualityScores) {
		covBiasType = coverageBiasType;
		clsfyType = classificationType;
		profileType = errorProfileType;

		if (profileType == ErrorProfileType::MACHINE_LEARNING) {
			ecu = ErrorCorrectionUnit(correctionType, epuClassify, kmerClassifier, correctIndels);
		} else if (profileType == ErrorProfileType::MOTIF_STATS_ONLY) {
			ecu = ErrorCorrectionUnit(correctionType, epuMotif, kmerClassifier, correctIndels);
		} else {
			ecu = ErrorCorrectionUnit(correctionType, epuOverall, kmerClassifier, correctIndels);
		}

	}

	void learnCoverageBias() {
		if (covBiasType == CoverageBiasType::IGNORE_BIAS) {
			return;
		}

		if (covBiasType == CoverageBiasType::MEDIAN_BIAS_REFERENCE) {
			biasTypeTitle = "median_bias_reference";
		} else if (covBiasType == CoverageBiasType::MEDIAN_BIAS_ALIGNMENT) {
			biasTypeTitle = "median_bias_alignment";
		} else if (covBiasType == CoverageBiasType::MEDIAN_BIAS_READS_ONLY) {
			biasTypeTitle = "median_bias_reads";
		}

		std::ifstream infile(dataset.plotPath + "coverageBias." + biasTypeTitle + ".txt");
		if (infile.good()) {
			std::cout << "Coverage bias has been already learned. Loading learned values instead.\n";
			covBias.loadBias(dataset.plotPath + "coverageBias." + biasTypeTitle + ".txt");
			return;
		}

		if (covBiasType == CoverageBiasType::MEDIAN_BIAS_REFERENCE) {
			covBias.learnBiasFromReferenceMatches(dataset.genome, counterReference, counterReads, dataset.readLengths);
		} else if (covBiasType == CoverageBiasType::MEDIAN_BIAS_ALIGNMENT) {
			covBias.learnBiasFromReferenceAlignment(dataset.genome, counterReference, dataset.readAlignmentsFileName,
					dataset.readLengths);
		} else if (covBiasType == CoverageBiasType::MEDIAN_BIAS_READS_ONLY) {
			covBias.learnBiasFromReadsOnly(dataset.readsOnlyFileName, counterReads, dataset.readLengths);
		}
		covBias.storeBias(dataset.plotPath + "coverageBias." + biasTypeTitle + ".txt");
		covBias.printBias();
		covBias.plotBias(dataset.plotPath + "coverageBiasPlot." + biasTypeTitle);
	}

	void trainKmerClassification() {
		if (clsfyType != KmerClassificationType::CLASSIFICATION_MACHINE_LEARNING) {
			return;
		} else {
			std::ifstream infile(dataset.plotPath + "bestKmerClassifier.joblib.pkl");
			if (infile.good()) {
				std::cout
						<< "Machine Learning K-mer Classification has been already learned. Loading learned values instead.\n";
				kmerClassifier.loadClassifier(dataset.plotPath + "bestKmerClassifier.joblib.pkl");
			} else {
				std::cout << "Training machine learning - Kmer-Classifier...\n";
				kmerClassifier.trainClassifier(dataset, counterReference);
				std::cout << "Finished training machine learning - Kmer-Classifier.\n";
				std::cout << "Storing machine learning - Kmer-Classifier...\n";
				kmerClassifier.storeClassifier(dataset.plotPath + "bestKmerClassifier.joblib.pkl");
				std::cout << "Finished storing machine learning - Kmer-Classifier.\n";
			}
		}
	}

	void extractErrors() {
		std::ifstream infile(dataset.readAlignmentsFileName + ".trueCorrections.txt");
		if (infile.good()) {
			std::cout << "Errors have already been extracted. Skipping error extraction.\n";
		} else {
			edu.addAlignmentsFile(dataset.readAlignmentsFileName, dataset.plotPath);
			edu.correctReads(dataset.genome);
			//edu.correctReadsMultithreaded(dataset.genome);
		}
	}

	void trainErrorProfile() {
		if (profileType == ErrorProfileType::OVERALL_STATS_ONLY) {
			std::ifstream infile(dataset.plotPath + "overallErrorProfile.txt");
			if (infile.good()) {
				std::cout << "Overall Error Profile has been already learned. Loading learned values instead.\n";
				epuOverall.loadErrorProfile(dataset.plotPath + "overallErrorProfile.txt", counterReads);
			} else {
				std::cout << "Training overall error profile...\n";
				epuOverall.learnErrorProfileFromFilesAligned(dataset.readAlignmentsFileName + ".trueCorrections.txt");
				//epuOverall.plotErrorProfile();
				std::cout << "Finished training overall error profile.\n";
				std::cout << "Storing overall error profile...\n";
				epuOverall.storeErrorProfile(dataset.plotPath + "overallErrorProfile.txt");
				std::cout << "Finished storing overall error profile...\n";
			}
		} else if (profileType == ErrorProfileType::MOTIF_STATS_ONLY) {
			std::ifstream infile(dataset.plotPath + "motifErrorProfile.txt");
			if (infile.good()) {
				std::cout << "Motif Error Profile has been already learned. Loading learned values instead.\n";
				epuMotif.loadErrorProfile(dataset.plotPath + "motifErrorProfile.txt", counterReads);
			} else {
				std::cout << "Training motif error profile...\n";
				epuMotif.learnErrorProfileFromFilesAligned(dataset.readAlignmentsFileName + ".trueCorrections.txt");
				//epuMotif.plotErrorProfile();
				std::cout << "Finished training motif error profile.\n";
				std::cout << "Storing motif error profile...\n";
				epuMotif.storeErrorProfile(dataset.plotPath + "motifErrorProfile.txt");
				std::cout << "Finished storing motif error profile...\n";
			}
		} else { // the full machine learning machinery
			std::ifstream infile(dataset.plotPath + "motifErrorProfile.txt");
			if (infile.good()) {
				std::cout << "Motif Error Profile has been already learned. Loading learned values instead.\n";
				epuMotif.loadErrorProfile(dataset.plotPath + "motifErrorProfile.txt", counterReads);
			} else {
				std::cout << "Training motif error profile...\n";
				epuMotif.learnErrorProfileFromFilesAligned(dataset.readAlignmentsFileName + ".trueCorrections.txt");
				//epuMotif.plotErrorProfile();
				std::cout << "Finished training motif error profile.\n";
				std::cout << "Storing motif error profile...\n";
				epuMotif.storeErrorProfile(dataset.plotPath + "motifErrorProfile.txt");
				std::cout << "Finished storing motif error profile...\n";
			}

			std::ifstream infile2(dataset.plotPath + biasTypeTitle + ".classifierErrorProfile.txt");
			if (infile2.good()) {
				std::cout << "Classifier Error Profile has been already learned. Loading learned values instead.\n";
				epuClassify.loadErrorProfile(dataset.plotPath + biasTypeTitle + ".classifierErrorProfile.txt", counterReads);
			} else {
				std::cout << "Training classifier error profile...\n";
				epuClassify.learnErrorProfileFromFilesAligned(dataset.readAlignmentsFileName + ".trueCorrections.txt", dataset.acceptProb);
				std::cout << "Finished training classifier error profile.\n";
				//epuClassify.plotErrorProfile();
				std::cout << "Storing classifier error profile...\n";
				epuClassify.storeErrorProfile(dataset.plotPath + biasTypeTitle + ".classifierErrorProfile.txt");
				std::cout << "Finished storing classifier error profile.\n";
			}
		}
	}

	void correctReads() {
		ecu.addReadsFile(dataset.readsFileName, dataset.plotPath);
		std::cout << "Correcting reads...\n";
		//ecu.correctReadsMultithreaded();
		ecu.correctReads();
		std::cout << "Finished correcting reads.\n";
	}

	Dataset dataset;
	KmerCounter counterReads;
	KmerCounter counterReference;
	PerfectUniformSequencingModel pusm;

	CoverageBiasUnit covBias;

	KmerClassificationUnit kmerClassifier;

	OverallErrorProfile epuOverall;
	MotifErrorProfile epuMotif;
	ClassifierErrorProfile epuClassify;

	ErrorDetectionUnit edu;
	ErrorCorrectionUnit ecu;
private:
	CoverageBiasType covBiasType;
	KmerClassificationType clsfyType;
	ErrorProfileType profileType;

	std::string biasTypeTitle;
};
