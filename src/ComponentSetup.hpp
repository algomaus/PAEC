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
#include "ErrorCorrectionEvaluation.h"

#include "external/gnuplot-iostream.h"

using namespace std::placeholders;

enum ErrorProfileType {
	OVERALL_STATS_ONLY = 0, MOTIF_STATS_ONLY = 1, MACHINE_LEARNING = 2
};

class ComponentSetup {
public:
	ComponentSetup(Dataset &ds, CoverageBiasType coverageBiasType, KmerClassificationType classificationType,
			ErrorProfileType errorProfileType, ErrorCorrectionType corrType, bool indels) :
			dataset(ds), counterReads(ds.readsOnlyFileName), counterReference(ds.referenceFileName), pusm(
					ds.readLengthsPtr, ds.genomeSize, ds.genomeType), covBias(coverageBiasType, ds.genomeSize, pusm), kmerClassifier(
					counterReads, counterReference, covBias, pusm, classificationType), epuMotif(counterReads), epuMotif2(counterReads), epuClassify(
					ds.plotPath, coverageBiasType, kmerClassifier, epuMotif, ds.hasQualityScores), epuClassify2(
					ds.plotPath, coverageBiasType, kmerClassifier, epuMotif2, ds.hasQualityScores) {
		covBiasType = coverageBiasType;
		clsfyType = classificationType;
		profileType = errorProfileType;
		correctionType = corrType;
		correctIndels = indels;

		edu = ErrorDetectionUnit(ece);

		if (profileType == ErrorProfileType::MACHINE_LEARNING) {
			ecu = ErrorCorrectionUnit(correctionType, epuClassify, kmerClassifier, correctIndels, ece);
			//edu.addObserver(epuMotif);
			//edu.addObserver(epuClassify);
			//ecu.addObserver(epuMotif2);
			//ecu.addObserver(epuClassify2);
			ecu.addObserver(epuOverall2);
		} else if (profileType == ErrorProfileType::MOTIF_STATS_ONLY) {
			ecu = ErrorCorrectionUnit(correctionType, epuMotif, kmerClassifier, correctIndels, ece);
			//edu.addObserver(epuMotif);
			//ecu.addObserver(epuMotif2);
		} else {
			ecu = ErrorCorrectionUnit(correctionType, epuOverall, kmerClassifier, correctIndels,ece);
			//edu.addObserver(epuOverall);
			//ecu.addObserver(epuOverall2);
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
			covBias.learnBiasFromReadsOnly(dataset.readsFileName, counterReads, dataset.readLengths);
		}
		covBias.storeBias(dataset.plotPath + "coverageBias." + biasTypeTitle + ".txt");
		covBias.printBias();
		covBias.plotBias(dataset.plotPath + "coverageBiasPlot." + biasTypeTitle);
	}

	void experimentAllErrorProfiles() {
		CoverageBiasUnit biasUnitLoaded(CoverageBiasType::MEDIAN_BIAS_READS_ONLY, dataset.genomeSize, pusm);
		KmerClassificationUnit classifier(counterReads, counterReference, biasUnitLoaded, pusm,
				KmerClassificationType::CLASSIFICATION_MACHINE_LEARNING);
		classifier.loadClassifier(dataset.plotPath + "bestKmerClassifier.joblib.pkl");

		// infer overall error profile
		OverallErrorProfile overallError;

		// extract training data
		ErrorDetectionUnit detect;
		detect.addAlignmentsFile(dataset.readAlignmentsFileName, dataset.plotPath);
		detect.addObserver(overallError);
		detect.correctReads(dataset.genome);
		//overallError.learnErrorProfileFromFilesAligned(dataset.plotPath + "trueCorrections.txt");
		overallError.plotErrorProfile();
		std::cout << "Finished training overall error profile.\n";
		std::cout << "Storing overall error profile...\n";
		overallError.storeErrorProfile(dataset.plotPath + "overallErrorProfile.txt");
		std::cout << "Finished storing overall error profile...\n";

		// infer motif error profile
		MotifErrorProfile motifError(counterReads);

		std::cout << "Training motif error profile...\n";
		ErrorDetectionUnit detect2;
		detect2.addAlignmentsFile(dataset.readAlignmentsFileName, dataset.plotPath);
		detect2.addObserver(motifError);
		detect2.correctReads(dataset.genome);
		//motifError.learnErrorProfileFromFilesAligned(dataset.plotPath + "trueCorrections.txt");
		motifError.plotErrorProfile();
		std::cout << "Finished training motif error profile.\n";
		std::cout << "Storing motif error profile...\n";
		motifError.storeErrorProfile(dataset.plotPath + "motifErrorProfile.txt");
		std::cout << "Finished storing motif error profile...\n";

		ClassifierErrorProfile classi(dataset.plotPath, CoverageBiasType::MEDIAN_BIAS_READS_ONLY, classifier,
				motifError, true);
		// infer machine learning error profile
		std::cout << "Training classifier error profile...\n";
		ErrorDetectionUnit detect3;
		detect3.addAlignmentsFile(dataset.readAlignmentsFileName, dataset.plotPath);
		detect3.addObserver(classi);
		detect3.correctReads(dataset.genome);
		//classi.learnErrorProfileFromFilesAligned(dataset.plotPath + "trueCorrections.txt", dataset.acceptProb);
		std::cout << "Finished training classifier error profile.\n";
		//classi.plotErrorProfile();
		std::cout << "Storing classifier error profile...\n";
		classi.storeErrorProfile(dataset.plotPath + biasTypeTitle + ".classifierErrorProfile.txt");
		std::cout << "Finished storing classifier error profile.\n";
	}

	void experimentAllKmerClassifiers() {
		CoverageBiasUnit biasUnitLoaded(CoverageBiasType::MEDIAN_BIAS_READS_ONLY, dataset.genomeSize, pusm);
		biasUnitLoaded.loadBias(dataset.plotPath + "coverageBias.median_bias_reads.txt");
		KmerClassificationUnit classifier(counterReads, counterReference, biasUnitLoaded, pusm,
				KmerClassificationType::CLASSIFICATION_MACHINE_LEARNING);
		classifier.trainClassifier(dataset, counterReference);
		classifier.storeClassifier(dataset.plotPath + "bestKmerClassifier.joblib.pkl");
	}

	void experimentAllCoverageBiases() {
		std::string biasTitle = "median_bias_reference";
		CoverageBiasUnit covBiasAct(CoverageBiasType::MEDIAN_BIAS_REFERENCE, dataset.genomeSize, pusm);
		covBiasAct.learnBiasFromReferenceMatches(dataset.genome, counterReference, counterReads, dataset.readLengths);
		covBiasAct.storeBias(dataset.plotPath + "coverageBias." + biasTitle + ".txt");
		covBiasAct.printBias();
		covBiasAct.plotBias(dataset.plotPath + "coverageBiasPlot." + biasTitle);

		biasTitle = "median_bias_alignment";
		CoverageBiasUnit covBiasAct2(CoverageBiasType::MEDIAN_BIAS_ALIGNMENT, dataset.genomeSize, pusm);
		covBiasAct2.learnBiasFromReferenceAlignment(dataset.genome, counterReference, dataset.readAlignmentsFileName,
				dataset.readLengths);
		covBiasAct2.storeBias(dataset.plotPath + "coverageBias." + biasTitle + ".txt");
		covBiasAct2.printBias();
		covBiasAct2.plotBias(dataset.plotPath + "coverageBiasPlot." + biasTitle);

		biasTitle = "median_bias_reads";
		CoverageBiasUnit covBiasAct3(CoverageBiasType::MEDIAN_BIAS_READS_ONLY, dataset.genomeSize, pusm);
		covBiasAct3.learnBiasFromReadsOnly(dataset.readsFileName, counterReads, dataset.readLengths);
		covBiasAct3.storeBias(dataset.plotPath + "coverageBias." + biasTitle + ".txt");
		covBiasAct3.printBias();
		covBiasAct3.plotBias(dataset.plotPath + "coverageBiasPlot." + biasTitle);

		std::string filename = dataset.plotPath + "medianCoverageBiases";
		Gnuplot gp("tee " + filename + ".gnu | gnuplot -persist");
		double minX = 0.0;
		double maxX = 1.0;
		double minY = 0.0;
		double maxY = 2.0;

		gp << "set terminal png size 800,600 enhanced font \"Helvetica,10\"" << std::endl;
		gp << "set output '" << filename << ".png'" << std::endl;
		gp << "set xrange [" << minX << ":" << maxX << "]\nset yrange [" << minY << ":" << maxY << "]\n";
		gp << "set xlabel \"" << "G/C content" << "\"" << std::endl;
		gp << "set ylabel \"" << "coverage bias" << "\"" << std::endl;

		std::string titlePart = dataset.name;
		for (size_t i = 0; i < titlePart.size(); ++i) {
			if (titlePart[i] == '_') {
				titlePart[i] = ' ';
			}
		}

		gp << "set title \"" << titlePart << " Median Coverage Biases\"" << std::endl;
		gp << "set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 1.5" << std::endl;
		gp << "set style line 2 lc rgb '#00994c' lt 1 lw 2 pt 7 ps 1.5" << std::endl;
		gp << "set style line 3 lc rgb '#99004c' lt 1 lw 2 pt 7 ps 1.5" << std::endl;
		gp << "set datafile separator \",\"" << std::endl;

		gp << "plot '" << dataset.plotPath
				<< "coverageBiasPlot.median_bias_reads_data.csv' with linespoints ls 1 title 'reads only',\\\n'"
				<< dataset.plotPath
				<< "coverageBiasPlot.median_bias_alignment_data.csv' with linespoints ls 2 title 'alignment',\\\n'"
				<< dataset.plotPath
				<< "coverageBiasPlot.median_bias_reference_data.csv' with linespoints ls 3 title 'reference'"
				<< std::endl;
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
		std::ifstream infile(dataset.plotPath + "trueCorrections.txt");
		if (infile.good()) {
			std::cout << "Errors have already been extracted. Skipping error extraction.\n";
		} else {
			edu.addAlignmentsFile(dataset.readAlignmentsFileName, dataset.plotPath);
			edu.correctReads(dataset.genome);
			//edu.correctReadsMultithreaded(dataset.genome);

			if (profileType == ErrorProfileType::OVERALL_STATS_ONLY) {
				epuOverall.plotErrorProfile();
			}
		}
	}

	void trainErrorProfile() {
		if (profileType == ErrorProfileType::OVERALL_STATS_ONLY) {
			std::ifstream infile(dataset.plotPath + "overallErrorProfile.txt");
			if (infile.good()) {
				std::cout << "Overall Error Profile has been already learned. Loading learned values instead.\n";
				epuOverall.loadErrorProfile(dataset.plotPath + "overallErrorProfile.txt", counterReads);
				epuOverall.plotErrorProfile();
			} else {
				std::cout << "Training overall error profile...\n";
				epuOverall.learnErrorProfileFromFilesAligned(dataset.plotPath + "trueCorrections.txt");
				epuOverall.plotErrorProfile();
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
				epuMotif.learnErrorProfileFromFilesAligned(dataset.plotPath + "trueCorrections.txt");
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
				epuMotif.learnErrorProfileFromFilesAligned(dataset.plotPath + "trueCorrections.txt");
				//epuMotif.plotErrorProfile();
				std::cout << "Finished training motif error profile.\n";
				std::cout << "Storing motif error profile...\n";
				epuMotif.storeErrorProfile(dataset.plotPath + "motifErrorProfile.txt");
				std::cout << "Finished storing motif error profile...\n";
			}

			std::ifstream infile2(dataset.plotPath + ".classifierErrorProfile.txt");
			if (infile2.good()) {
				std::cout << "Classifier Error Profile has been already learned. Loading learned values instead.\n";
				epuClassify.loadErrorProfile(dataset.plotPath + ".classifierErrorProfile.txt",
						counterReads);
			} else {
				std::cout << "Could not find " << dataset.plotPath + ".classifierErrorProfile.txt" << "\n";

				std::cout << "Training classifier error profile...\n";
				epuClassify.learnErrorProfileFromFilesAligned(dataset.plotPath + "trueCorrections.txt",
						dataset.acceptProb);
				std::cout << "Finished training classifier error profile.\n";
				//epuClassify.plotErrorProfile();
				std::cout << "Storing classifier error profile...\n";
				epuClassify.storeErrorProfile(dataset.plotPath + ".classifierErrorProfile.txt");
				std::cout << "Finished storing classifier error profile.\n";
			}
		}
	}

	void correctReads() {
		ecu.addReadsFile(dataset.readsFileName, dataset.plotPath);
		std::cout << "Correcting reads, Part 1...\n";
		//ecu.correctReadsMultithreaded();
		ecu.correctReads();
		std::cout << "Finished correcting reads.\n";

		/*
		std::cout << "Correcting reads, Part 2...\n";
		if (profileType == ErrorProfileType::MACHINE_LEARNING) {
			ecu = ErrorCorrectionUnit(correctionType, epuClassify2, kmerClassifier, correctIndels);
		} else if (profileType == ErrorProfileType::MOTIF_STATS_ONLY) {
			ecu = ErrorCorrectionUnit(correctionType, epuMotif2, kmerClassifier, correctIndels);
		} else {
			epuOverall2.plotErrorProfile();
			ecu = ErrorCorrectionUnit(correctionType, epuOverall2, kmerClassifier, correctIndels);
		}
		ecu.addReadsFile(dataset.readsFileName, dataset.plotPath);
		//ecu.correctReadsMultithreaded();
		ecu.correctReads();
		std::cout << "Finished correcting reads, Part 2.\n";
		*/
	}

	Dataset dataset;
	KmerCounter counterReads;
	KmerCounter counterReference;
	PerfectUniformSequencingModel pusm;

	CoverageBiasUnit covBias;

	KmerClassificationUnit kmerClassifier;

	OverallErrorProfile epuOverall, epuOverall2;
	MotifErrorProfile epuMotif, epuMotif2;
	ClassifierErrorProfile epuClassify, epuClassify2;

	ErrorDetectionUnit edu;
	ErrorCorrectionUnit ecu;
	ErrorCorrectionEvaluation ece;
private:
	CoverageBiasType covBiasType;
	KmerClassificationType clsfyType;
	ErrorProfileType profileType;
	ErrorCorrectionType correctionType;
	bool correctIndels;

	std::string biasTypeTitle;
};
