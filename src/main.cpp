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

void initPython() {
	Py_Initialize();

	PyObject *sys = PyImport_ImportModule("sys");
	PyObject *path = PyObject_GetAttrString(sys, "path");
	PyList_Append(path,
			PyString_FromString("/home/sarah/Documents/Master Thesis Topic Extension/PAEC/src/python"));
	Py_DECREF(sys);
	Py_DECREF(path);
}

void exitPython() {
	Py_Finalize();
}

void processDataset(Dataset &ds) {
	bool extendedCover = false; // do Step 2 of the error correction method or not
	ComponentSetup cs(ds, CoverageBiasType::MEDIAN_BIAS_READS_ONLY,
			KmerClassificationType::CLASSIFICATION_CHEATING, ErrorProfileType::MACHINE_LEARNING,
			ErrorCorrectionType::KMER_BASED, extendedCover);

	// extract errors from alignment
	cs.extractErrors();

	// learn coverage bias, or load it if it is already there.
	cs.learnCoverageBias();

	// learn coverage bias.
	//cs.learnCoverageBias();
	// train machine learning method for k-mer classification
	//cs.trainKmerClassification();

	// train error profile
	cs.trainErrorProfile();

	// correct reads
	cs.correctReads();

	/*
	// train machine learning method for k-mer classification
	cs.trainKmerClassification();



	*/

	// *** THE EXPERIMENTS ***
	//cs.experimentAllCoverageBiases();
	//cs.experimentAllKmerClassifiers();
	//cs.experimentAllErrorProfiles();
}

int main() {
	initPython();

	std::vector<Dataset> datasets;
	//datasets.push_back(Dataset("ecoli_illumina_simulated"));
	//datasets.push_back(Dataset("SRR396536"));
	datasets.push_back(Dataset("ebola_illumina_simulated"));
	//datasets.push_back(Dataset("ecoli_pacbio_simulated"));
	//datasets.push_back(Dataset("SRR396537"));
	//datasets.push_back(Dataset("ebola_pacbio_simulated"));
	//datasets.push_back(Dataset("SRR1284073"));

	for (size_t i = 0; i < datasets.size(); ++i) {
		std::cout << "\n" << datasets[i].name << ":\n";
		processDataset(datasets[i]);
	}

	exitPython();
	return 0;
}
