//============================================================================
// Name        : main.cpp
// Author      : Sarah Lutteropp
// Version     :
// Copyright   : Your copyright notice
// Description :
//============================================================================

#include <chrono>
#include <iostream>
#include <memory>
#include <string>
#include <python2.7/Python.h>

#include "AlignedInformation/Dataset.hpp"
#include "AlignedInformation/ErrorDetectionUnit.h"
#include "AlignedInformation/GenomeReader.h"
#include "ComponentSetup.hpp"
#include "CoverageBias/CoverageBiasUnit.h"
#include "CoverageBias/PUSM.h"
#include "ErrorProfile/ErrorProfileUnit.hpp"
#include "ErrorProfile/machine_learning/ClassifierErrorProfile.h"
#include "ErrorProfile/motif_analysis/MotifErrorProfile.h"
#include "ErrorProfile/OverallErrorProfile.h"
#include "KmerClassification/KmerClassificationUnit.h"
#include "KmerClassification/KmerCounter.h"

using namespace std;

/*void detectErrorsFromAlignments(Dataset &ds) {
 ErrorDetectionUnit edu;

 edu.addAlignmentsFile(ds.readAlignmentsFileName, ds.plotPath);

 std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
 edu.correctReadsMultithreaded(ds.genome);
 //edu.correctReads(ds.genome);

 std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
 std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()
 << std::endl;
 }

 void learnManyProfilesFromAlignments(Dataset &ds, bool useCoverageBias) {
 std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

 std::shared_ptr<KmerCounter> counter = std::make_shared<KmerCounter>(ds.readsFileName);
 OverallErrorProfile oep;
 oep.learnErrorProfileFromFilesAligned(ds.plotPath + "corrections.txt");
 oep.storeErrorProfile(ds.plotPath + "overallProfile.txt");
 oep.plotErrorProfile();

 std::shared_ptr<MotifErrorProfile> mep = std::make_shared<MotifErrorProfile>(counter);
 mep->learnErrorProfileFromFilesAligned(ds.plotPath + "corrections.txt");
 mep->storeErrorProfile(ds.plotPath + "motifProfile.txt");

 std::shared_ptr<PerfectUniformSequencingModel> pusm = std::make_shared<PerfectUniformSequencingModel>(
 ds.readLengthsPtr, ds.genomeSize, ds.genomeType);

 std::shared_ptr<CoverageBiasUnit> covBias;
 std::shared_ptr<KmerClassificationUnit> kmerClassifier;
 if (!useCoverageBias) {
 covBias = std::make_shared<CoverageBiasUnit>(CoverageBiasType::IGNORE_BIAS, ds.genomeSize, pusm); // no coverage bias
 kmerClassifier = std::make_shared<KmerClassificationUnit>(counter, covBias, pusm,
 KmerClassificationType::CLASSIFICATION_STATISTICAL);
 } else {
 covBias = std::make_shared<CoverageBiasUnit>(CoverageBiasType::MEDIAN_BIAS_REFERENCE, ds.genomeSize, pusm);
 std::shared_ptr<KmerCounter> counterRef = std::make_shared<KmerCounter>(ds.referenceFileName);
 std::string referenceGenome = GenomeReader::readGenome(ds.referenceFileName);
 covBias->learnBiasFromReferenceMatches(ds.genome, counterRef, counter, ds.readLengths);
 //covBias->learnBiasFromReference(ds.genome, counterRef, ds.readAlignmentsFileName, ds.readLengths);
 covBias->storeBias(ds.plotPath + "coverageBias.txt");
 covBias->plotBias();
 kmerClassifier = std::make_shared<KmerClassificationUnit>(counter, covBias, pusm,
 KmerClassificationType::CLASSIFICATION_STATISTICAL);

 }

 ClassifierErrorProfile cep(ds.plotPath, kmerClassifier, mep, true);
 cep.learnErrorProfileFromFilesAligned(ds.plotPath + "corrections.txt");
 cep.storeErrorProfile(ds.plotPath + "mlearnProfile.txt");

 std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
 std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()
 << std::endl;
 }*/

void initPython() {
	Py_Initialize();

	PyObject *sys = PyImport_ImportModule("sys");
	PyObject *path = PyObject_GetAttrString(sys, "path");
	PyList_Append(path,
			PyString_FromString("/home/sarah/Documents/Masterarbeit/devel/MasterThesisCodeV2/PAEC/src/python"));
	Py_DECREF(sys);
	Py_DECREF(path);
}

void exitPython() {
	Py_Finalize();
}

void processDataset(Dataset &ds) {
	// TODO: FIXME: Bug in predict proba with classifier profile
	ComponentSetup cs(ds, CoverageBiasType::MEDIAN_BIAS_REFERENCE,
			KmerClassificationType::CLASSIFICATION_MACHINE_LEARNING, ErrorProfileType::MACHINE_LEARNING,
			ErrorCorrectionType::KMER_BASED);
	// learn coverage bias.
	cs.learnCoverageBias();
	// train machine learning method for k-mer classification
	cs.trainKmerClassification();

	// extract errors from alignment
	cs.extractErrors();
	// train error profile
	cs.trainErrorProfile();
	// precorrect reads
	cs.precorrectReads();
	// postcorrect reads
	cs.postcorrectReads();
}

int main() {
	initPython();

	Dataset ds("ebola_simulated_pacbio");
	//Dataset ds("ebola_simulated_illumina");
	//Dataset ds("SRR396537");

	//detectErrorsFromAlignments(ds);
	//learnManyProfilesFromAlignments(ds, true);
	processDataset(ds);

	exitPython();
	return 0;
}
