/*
 * ClassifierErrorProfile.cpp
 *
 *  Created on: Feb 8, 2017
 *      Author: sarah
 */

#include "ClassifierErrorProfile.h"
#include <cassert>
#include <fstream>
#include <random>
#include "../../AlignedInformation/CorrectionAligned.h"

ClassifierErrorProfile::ClassifierErrorProfile() {
	finalized = false;
	classifierNextGap = NULL;
	classifierCurrentBase = NULL;
	overallErrorRateCurrentBase = 0.01;
	overallErrorRateNextGap = 0.01;
	useQual = true;
	alreadyHasTrainingData = false;
	useKmerZScores = false;
	errorsOnly = true;
}

void ClassifierErrorProfile::setOverallErrorRates(double errorRateBase, double errorRateGap) {
	overallErrorRateCurrentBase = errorRateBase;
	overallErrorRateNextGap = errorRateGap;
}

ClassifierErrorProfile::ClassifierErrorProfile(const std::string &plotPath, CoverageBiasType coverageBiasType,
		KmerClassificationUnit &kmerClassifier, MotifErrorProfile &motifProfile, bool useQualityScores,
		bool useKmerZScoresOrNot, bool useErrorsOnly) {
	finalized = false;
	useQual = useQualityScores;
	useKmerZScores = useKmerZScoresOrNot;
	errorsOnly = useErrorsOnly;

	std::string biasTypeTitle;
	if (coverageBiasType == CoverageBiasType::MEDIAN_BIAS_REFERENCE) {
		biasTypeTitle = "median_bias_reference";
	} else if (coverageBiasType == CoverageBiasType::MEDIAN_BIAS_ALIGNMENT) {
		biasTypeTitle = "median_bias_alignment";
	} else if (coverageBiasType == CoverageBiasType::MEDIAN_BIAS_READS_ONLY) {
		biasTypeTitle = "median_bias_reads";
	}

	clsfyCurrentBase = plotPath + biasTypeTitle + ".clsfyCurrentBase.txt";
	trainCurrentBase = plotPath + biasTypeTitle + ".trainCurrentBase.csv";
	clsfyNextGap = plotPath + biasTypeTitle + ".clsfyNextGap.txt";
	trainNextGap = plotPath + biasTypeTitle + ".trainNextGap.csv";
	overallErrorRateCurrentBase = 0.1;
	overallErrorRateNextGap = 0.1;
	feCurrentBase = std::make_unique<FeatureExtractorCurrentBase>(kmerClassifier, motifProfile, useKmerZScores,
			errorsOnly);
	feNextGap = std::make_unique<FeatureExtractorNextGap>(kmerClassifier, motifProfile, useKmerZScores, errorsOnly);

	alreadyHasTrainingData = false;
	std::ifstream test1(trainCurrentBase);
	std::ifstream test2(trainNextGap);
	if (test1.good() && test2.good()) {
		alreadyHasTrainingData = true;
		std::cout << "Training data files for the ClassifierErrorProfile are already there! Great, skipping this part then, if needed.\n";
	}

	if (!alreadyHasTrainingData) {
		std::ofstream outfileCurrent(trainCurrentBase);
		if (!outfileCurrent.good()) {
			throw std::runtime_error("Could not create file: " + trainCurrentBase);
		}

		if (useQual) {
			outfileCurrent << feCurrentBase->getFeatureNames() << "\n";
		} else {
			outfileCurrent << feCurrentBase->getFeatureNamesNoQual() << "\n";
		}
		outfileCurrent.close();
		std::ofstream outfileNext(trainNextGap);
		if (!outfileNext.good()) {
			throw std::runtime_error("Could not create file: " + trainNextGap);
		}
		if (useQual) {
			outfileNext << feNextGap->getFeatureNames() << "\n";
		} else {
			outfileNext << feNextGap->getFeatureNamesNoQual() << "\n";
		}
		outfileNext.close();
	}

	distribution = std::uniform_real_distribution<double>(0.0, 1.0);

	/*Py_Initialize();

	 PyObject *sys = PyImport_ImportModule("sys");
	 PyObject *path = PyObject_GetAttrString(sys, "path");
	 PyList_Append(path,
	 PyString_FromString("/home/sarah/Documents/Master Thesis Topic Extension/PAEC/src/python"));
	 Py_DECREF(sys);
	 Py_DECREF(path);*/

	PyObject* module = PyImport_ImportModule((char*) "blackbox");
	assert(module != NULL);
	PyObject* klass = PyObject_GetAttrString(module, (char*) "classifier");
	assert(klass != NULL);
	classifierNextGap = PyInstance_New(klass, NULL, NULL);
	assert(classifierNextGap != NULL);
	classifierCurrentBase = PyInstance_New(klass, NULL, NULL);
	assert(classifierCurrentBase != NULL);

	Py_DECREF(klass);
	Py_DECREF(module);

	// set the feature names
	std::vector<std::string> featureNamesNextGap;
	std::vector<std::string> featureNamesCurrentBase;
	if (useQual) {
		featureNamesNextGap = feNextGap->getFeatureNamesVector();
		featureNamesCurrentBase = feCurrentBase->getFeatureNamesVector();
	} else {
		featureNamesNextGap = feNextGap->getFeatureNamesNoQualVector();
		featureNamesCurrentBase = feCurrentBase->getFeatureNamesNoQualVector();
	}
	// convert the feature name vectors into python objects
	PyObject* pyFeatureNamesNextGap = PyList_New(featureNamesNextGap.size());
	for (size_t i = 0; i < featureNamesNextGap.size(); ++i) {
		PyList_SetItem(pyFeatureNamesNextGap, i, PyString_FromString(featureNamesNextGap[i].c_str()));
	}
	PyObject* pyFeatureNamesCurrentBase = PyList_New(featureNamesCurrentBase.size());
	for (size_t i = 0; i < featureNamesCurrentBase.size(); ++i) {
		PyList_SetItem(pyFeatureNamesCurrentBase, i, PyString_FromString(featureNamesCurrentBase[i].c_str()));
	}
	PyObject* pyResultCurrent = PyObject_CallMethod(classifierNextGap, (char*) "set_features", (char*) "O",
			pyFeatureNamesNextGap);
	if (!pyResultCurrent) {
		throw std::runtime_error("set_features did not work");
	}
	PyObject* pyResultNext = PyObject_CallMethod(classifierCurrentBase, (char*) "set_features", (char*) "O",
			pyFeatureNamesCurrentBase);
	if (!pyResultNext) {
		throw std::runtime_error("set_features did not work");
	}

	// set the class id's
	PyObject* pyClassNamesNextGap = PyList_New(feNextGap->getClasses().size());
	PyObject* pyClassNamesCurrentBase = PyList_New(feCurrentBase->getClasses().size());
	// convert the class name vectors into python objects

	int idx = 0;
	for (size_t i = 0; i < feNextGap->getClasses().size(); ++i) {
		PyList_SetItem(pyClassNamesNextGap, idx, PyInt_FromLong(errorTypeToNumber(feNextGap->getClasses()[i])));
		idx++;
	}
	idx = 0;
	for (size_t i = 0; i < feCurrentBase->getClasses().size(); ++i) {
		PyList_SetItem(pyClassNamesCurrentBase, idx, PyInt_FromLong(errorTypeToNumber(feCurrentBase->getClasses()[i])));
		idx++;
	}
	pyResultCurrent = PyObject_CallMethod(classifierNextGap, (char*) "set_classes", (char*) "O", pyClassNamesNextGap);
	pyResultNext = PyObject_CallMethod(classifierCurrentBase, (char*) "set_classes", (char*) "O",
			pyClassNamesCurrentBase);
}

ClassifierErrorProfile::~ClassifierErrorProfile() {
	if (classifierCurrentBase != NULL) {
		Py_DECREF(classifierCurrentBase);
	}
	if (classifierNextGap != NULL) {
		Py_DECREF(classifierNextGap);
	}
	//Py_Finalize();
}

std::unordered_map<ErrorType, double> ClassifierErrorProfile::getKmerErrorProbabilities(const std::string &kmer,
		size_t positionInKmer) {
	if (kmer.find("_") != std::string::npos) {
		throw std::runtime_error("Invalid k-mer!");
	}
	if (useQual) {
		throw std::runtime_error("This is the wrong classifier for k-mer error probabilities!");
	}
	if (!finalized) {
		finalize();
	}
	return getErrorProbabilitiesFinalized(kmer, positionInKmer);
}

std::unordered_map<ErrorType, double> ClassifierErrorProfile::getErrorProbabilities(const FASTQRead &read,
		size_t positionInRead) {
	if (read.sequence[positionInRead] == '_') {
		throw std::runtime_error("Multidel!");
	}
	if (!useQual) {
		throw std::runtime_error("This is the wrong classifier for FASTQ error probabilities!");
	}
	if (!finalized) {
		finalize();
	}
	return getErrorProbabilitiesFinalized(read, positionInRead);
}

void ClassifierErrorProfile::loadErrorProfile(const std::string &filepath, KmerCounter &counter) {
	std::ifstream infile(filepath, std::ios::binary);
	if (!infile.good()) {
		throw std::runtime_error("The file " + filepath + " does not exist!");
	}
	cereal::BinaryInputArchive iarchive(infile);
	ClassifierErrorProfile cep;
	iarchive(cep);
	clsfyCurrentBase = cep.clsfyCurrentBase;
	trainCurrentBase = cep.trainCurrentBase;
	clsfyNextGap = cep.clsfyNextGap;
	trainNextGap = cep.trainNextGap;
	overallErrorRateCurrentBase = cep.overallErrorRateCurrentBase;
	overallErrorRateNextGap = cep.overallErrorRateNextGap;
	finalized = cep.finalized;
	useKmerZScores = cep.useKmerZScores;

	std::string bestClassifierCurrentBaseFilepath = filepath + ".bestClassifier.currentBase.joblib.pkl";
	std::string bestClassifierNextGapFilepath = filepath + ".bestClassifier.nextGap.joblib.pkl";
	PyObject_CallMethod(classifierCurrentBase, (char*) "load_classifier", (char*) "s",
			bestClassifierCurrentBaseFilepath.c_str());
	PyObject_CallMethod(classifierNextGap, (char*) "load_classifier", (char*) "s",
			bestClassifierNextGapFilepath.c_str());
}
void ClassifierErrorProfile::storeErrorProfile(const std::string &filepath) {
	std::ofstream outfile(filepath, std::ios::binary);
	cereal::BinaryOutputArchive oarchive(outfile);
	oarchive(*this);

	std::string bestClassifierCurrentBaseFilepath = filepath + ".bestClassifier.currentBase.joblib.pkl";
	std::string bestClassifierNextGapFilepath = filepath + ".bestClassifier.nextGap.joblib.pkl";
	PyObject_CallMethod(classifierCurrentBase, (char*) "store_classifier", (char*) "s",
			bestClassifierCurrentBaseFilepath.c_str());
	PyObject_CallMethod(classifierNextGap, (char*) "store_classifier", (char*) "s",
			bestClassifierNextGapFilepath.c_str());
}
void ClassifierErrorProfile::plotErrorProfile() {
	// TODO ...
}

void ClassifierErrorProfile::reset() {
	finalized = false;
	if (alreadyHasTrainingData)
		return;
	std::ofstream outfileCurrent(trainCurrentBase);
	if (outfileCurrent.good()) {
		if (useQual) {
			outfileCurrent << feCurrentBase->getFeatureNames() << "\n";
		} else {
			outfileCurrent << feCurrentBase->getFeatureNamesNoQual() << "\n";
		}
		outfileCurrent.close();
	}
	std::ofstream outfileNext(trainNextGap);
	if (outfileNext.good()) {
		outfileNext << "";
		if (useQual) {
			outfileNext << feNextGap->getFeatureNames() << "\n";
		} else {
			outfileNext << feNextGap->getFeatureNamesNoQual() << "\n";
		}
		outfileNext.close();
	}
}

void ClassifierErrorProfile::processCorrection(const FASTQRead &read, size_t posInRead, ErrorType type,
		std::vector<bool> &nodel, std::vector<bool> &correct) {
	if (alreadyHasTrainingData)
		return;

	std::ofstream outCurrent(trainCurrentBase, std::ofstream::app);
	std::ofstream outNext(trainNextGap, std::ofstream::app);
	if (!outCurrent.good() || !outNext.good()) {
		throw std::runtime_error("The filepaths are not good.");
	}

	std::string trainingString;
	if (type == ErrorType::INSERTION || type == ErrorType::SUB_FROM_A || type == ErrorType::SUB_FROM_C
			|| type == ErrorType::SUB_FROM_G || type == ErrorType::SUB_FROM_T || type == ErrorType::CORRECT) {
		if (useQual) {
			trainingString = feCurrentBase->getTrainingString(read, posInRead, type) + "\n";
		} else {
			trainingString = feCurrentBase->getTrainingStringNoQual(read.sequence, posInRead, type) + "\n";
		}
		outCurrent << trainingString;

	} else {
		if (useQual) {
			trainingString = feNextGap->getTrainingString(read, posInRead, type) + "\n";
		} else {
			trainingString = feNextGap->getTrainingStringNoQual(read.sequence, posInRead, type) + "\n";
		}
		outNext << trainingString;
	}

	if (type == ErrorType::DEL_OF_A || type == ErrorType::DEL_OF_C || type == ErrorType::DEL_OF_G
			|| type == ErrorType::DEL_OF_T || type == ErrorType::MULTIDEL) {
		nodel[posInRead] = false;
	} else {
		assert(type != ErrorType::CORRECT && type != ErrorType::NODEL);
		correct[posInRead] = false;
	}
}

void ClassifierErrorProfile::processCorrect(const FASTQRead &read, size_t posInRead) {
	if (alreadyHasTrainingData || errorsOnly)
		return;

	std::ofstream outCurrent(trainCurrentBase, std::ofstream::app);
	if (!outCurrent.good()) {
		throw std::runtime_error("The filepath is not good.");
	}
	std::string trainingString;
	if (useQual) {
		trainingString = feCurrentBase->getTrainingString(read, posInRead, ErrorType::CORRECT) + "\n";
	} else {
		trainingString = feCurrentBase->getTrainingStringNoQual(read.sequence, posInRead, ErrorType::CORRECT) + "\n";
	}
	outCurrent << trainingString;
}

void ClassifierErrorProfile::processNodel(const FASTQRead &read, size_t posInRead) {
	if (alreadyHasTrainingData || errorsOnly)
		return;

	std::ofstream outNext(trainNextGap, std::ofstream::app);
	if (!outNext.good()) {
		throw std::runtime_error("The filepath is not good.");
	}
	std::string trainingString;
	if (useQual) {
		trainingString = feNextGap->getTrainingString(read, posInRead, ErrorType::NODEL) + "\n";
	} else {
		trainingString = feNextGap->getTrainingStringNoQual(read.sequence, posInRead, ErrorType::NODEL) + "\n";
	}
	outNext << trainingString;
}

void ClassifierErrorProfile::check(const CorrectedRead &corrRead, double acceptProb) {
	finalized = false;
	if (alreadyHasTrainingData)
		return;

	std::vector<bool> nodel(corrRead.originalRead.sequence.size(), true);
	std::vector<bool> correct(corrRead.originalRead.sequence.size(), true);
	for (Correction corr : corrRead.corrections) {
		double roll = distribution(generator);
		if (roll <= acceptProb) {
			processCorrection(corrRead.originalRead, corr.originalReadPos, corr.type, nodel, correct);
		}
	}

	if (!errorsOnly) {
		for (size_t i = 0; i < corrRead.originalRead.sequence.size(); ++i) {
			if (correct[i]) {
				double roll = distribution(generator);
				if (roll <= overallErrorRateCurrentBase) {
					roll = distribution(generator);
					if (roll <= acceptProb) {
						processCorrect(corrRead.originalRead, i);
					}
				}
			}
			if (nodel[i]) {
				double roll = distribution(generator);
				if (roll <= overallErrorRateNextGap) {
					roll = distribution(generator);
					if (roll <= acceptProb) {
						processNodel(corrRead.originalRead, i);
					}
				}
			}
		}
	}
}
void ClassifierErrorProfile::checkAligned(const CorrectedReadAligned &corrRead, double acceptProb) {
	finalized = false;
	if (alreadyHasTrainingData)
		return;

	std::vector<bool> nodel(corrRead.originalRead.sequence.size(), true);
	std::vector<bool> correct(corrRead.originalRead.sequence.size(), true);

	for (size_t i = 0; i < corrRead.alignedCorrections.size(); ++i) {
		double roll = distribution(generator);
		if (roll <= acceptProb) {
			const CorrectionAligned &ca = corrRead.alignedCorrections[i];
			processCorrection(corrRead.originalRead, ca.correction.originalReadPos, ca.correction.type, nodel, correct);
		}
	}

	if (!errorsOnly) {
		for (size_t i = 0; i < corrRead.originalRead.sequence.size(); ++i) {
			if (correct[i]) {
				double roll = distribution(generator);
				if (roll <= overallErrorRateCurrentBase) {
					roll = distribution(generator);
					if (roll <= acceptProb) {
						processCorrect(corrRead.originalRead, i);
					}
				}
			}
			if (nodel[i]) {
				double roll = distribution(generator);
				if (roll <= overallErrorRateNextGap) {
					roll = distribution(generator);
					if (roll <= acceptProb) {
						processNodel(corrRead.originalRead, i);
					}
				}
			}
		}
	}
}

void ClassifierErrorProfile::finalize() {
	PyObject* res = PyObject_CallMethod(classifierCurrentBase, (char*) "set_csv_file", (char*) "s",
			trainCurrentBase.c_str());
	if (!res) {
		throw std::runtime_error("classifierCurrentBase set_csv_file went wrong");
	}
	res = PyObject_CallMethod(classifierNextGap, (char*) "set_csv_file", (char*) "s", trainNextGap.c_str());
	if (!res) {
		throw std::runtime_error("classifierNextGap set_csv_file went wrong");
	}
	finalized = true;
}

std::unordered_map<ErrorType, double> ClassifierErrorProfile::getErrorProbabilitiesFinalized(const std::string &kmer,
		size_t positionInKmer) {
	std::unordered_map<ErrorType, double> probas;

	std::vector<double> featuresCurrentBase;
	featuresCurrentBase = feCurrentBase->getFeatureVectorNoQual(kmer, positionInKmer);
	std::vector<double> featuresNextGap;
	featuresNextGap = feNextGap->getFeatureVectorNoQual(kmer, positionInKmer);
	// convert the vectors into python objects
	PyObject* pyFeaturesCurrent = PyList_New(featuresCurrentBase.size());
	for (size_t i = 0; i < featuresCurrentBase.size(); ++i) {
		PyList_SetItem(pyFeaturesCurrent, i, PyFloat_FromDouble(featuresCurrentBase[i]));
	}
	PyObject* pyFeaturesNext = PyList_New(featuresNextGap.size());
	for (size_t i = 0; i < featuresNextGap.size(); ++i) {
		PyList_SetItem(pyFeaturesNext, i, PyFloat_FromDouble(featuresNextGap[i]));
	}
	// call the classifiers
	std::vector<double> logProbaCurrent;
	std::vector<double> logProbaNext;
	PyObject* pyResultCurrent = PyObject_CallMethod(classifierCurrentBase, (char*) "proba", (char*) "O",
			pyFeaturesCurrent);
	PyObject* pyResultNext = PyObject_CallMethod(classifierNextGap, (char*) "proba", (char*) "O", pyFeaturesNext);
	for (int i = 0; i < PyList_Size(pyResultCurrent); ++i) {
		logProbaCurrent.push_back(PyFloat_AsDouble(PyList_GetItem(pyResultCurrent, i)));
	}
	for (int i = 0; i < PyList_Size(pyResultNext); ++i) {
		logProbaNext.push_back(PyFloat_AsDouble(PyList_GetItem(pyResultNext, i)));
	}
	Py_DECREF(pyResultCurrent);
	Py_DECREF(pyResultNext);
	Py_DECREF(pyFeaturesCurrent);
	Py_DECREF(pyFeaturesNext);

	const std::vector<ErrorType> &classesCurrent = feCurrentBase->getClasses();
	assert(classesCurrent.size() == logProbaCurrent.size());
	for (size_t i = 0; i < classesCurrent.size(); ++i) {
		probas[classesCurrent[i]] = logProbaCurrent[i];
	}
	const std::vector<ErrorType> &classesNext = feNextGap->getClasses();
	assert(classesNext.size() == logProbaNext.size());
	for (size_t i = 0; i < classesNext.size(); ++i) {
		probas[classesNext[i]] = logProbaNext[i];
	}
	return probas;

}

std::unordered_map<ErrorType, double> ClassifierErrorProfile::getErrorProbabilitiesFinalized(const FASTQRead &read,
		size_t positionInRead) {
	std::unordered_map<ErrorType, double> probas;

	std::vector<double> featuresCurrentBase;
	if (useQual) {
		featuresCurrentBase = feCurrentBase->getFeatureVector(read, positionInRead);
	} else {
		featuresCurrentBase = feCurrentBase->getFeatureVectorNoQual(read.sequence, positionInRead);
	}
	std::vector<double> featuresNextGap;
	if (useQual) {
		featuresNextGap = feNextGap->getFeatureVector(read, positionInRead);
	} else {
		featuresNextGap = feNextGap->getFeatureVectorNoQual(read.sequence, positionInRead);
	}
	// convert the vectors into python objects
	PyObject* pyFeaturesCurrent = PyList_New(featuresCurrentBase.size());
	for (size_t i = 0; i < featuresCurrentBase.size(); ++i) {
		PyList_SetItem(pyFeaturesCurrent, i, PyFloat_FromDouble(featuresCurrentBase[i]));
	}
	PyObject* pyFeaturesNext = PyList_New(featuresNextGap.size());
	for (size_t i = 0; i < featuresNextGap.size(); ++i) {
		PyList_SetItem(pyFeaturesNext, i, PyFloat_FromDouble(featuresNextGap[i]));
	}
	// call the classifiers
	std::vector<double> logProbaCurrent;
	std::vector<double> logProbaNext;
	PyObject* pyResultCurrent = PyObject_CallMethod(classifierCurrentBase, (char*) "proba", (char*) "O",
			pyFeaturesCurrent);
	PyObject* pyResultNext = PyObject_CallMethod(classifierNextGap, (char*) "proba", (char*) "O", pyFeaturesNext);
	for (int i = 0; i < PyList_Size(pyResultCurrent); ++i) {
		logProbaCurrent.push_back(PyFloat_AsDouble(PyList_GetItem(pyResultCurrent, i)));
	}
	for (int i = 0; i < PyList_Size(pyResultNext); ++i) {
		logProbaNext.push_back(PyFloat_AsDouble(PyList_GetItem(pyResultNext, i)));
	}
	Py_DECREF(pyResultCurrent);
	Py_DECREF(pyResultNext);
	Py_DECREF(pyFeaturesCurrent);
	Py_DECREF(pyFeaturesNext);

	const std::vector<ErrorType> &classesCurrent = feCurrentBase->getClasses();
	assert(classesCurrent.size() == logProbaCurrent.size());
	for (size_t i = 0; i < classesCurrent.size(); ++i) {
		probas[classesCurrent[i]] = logProbaCurrent[i];
	}
	const std::vector<ErrorType> &classesNext = feNextGap->getClasses();
	assert(classesNext.size() == logProbaNext.size());
	for (size_t i = 0; i < classesNext.size(); ++i) {
		probas[classesNext[i]] = logProbaNext[i];
	}
	return probas;
}
