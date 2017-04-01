/*
 * KmerClassificationUnit.cpp
 *
 *  Created on: Jan 31, 2017
 *      Author: sarah
 */

#include "KmerClassificationUnit.h"

#include <cmath>
#include <fstream>

#include "../CoverageBias/PUSM.h"

KmerClassificationUnit::KmerClassificationUnit(KmerCounter &kmerCounter, CoverageBiasUnit &biasUnitRef,
		PerfectUniformSequencingModel &pusmRef, KmerClassificationType type) :
		counter(kmerCounter), biasUnit(biasUnitRef), pusm(pusmRef) {
	classificationType = type;
	mlClassifier = NULL;
	if (type == KmerClassificationType::CLASSIFICATION_MACHINE_LEARNING) {
		PyObject* module = PyImport_ImportModule("blackbox");
		assert(module != NULL);
		PyObject* klass = PyObject_GetAttrString(module, "classifier");
		assert(klass != NULL);
		mlClassifier = PyInstance_New(klass, NULL, NULL);
		assert(mlClassifier != NULL);
		Py_DECREF(klass);
		Py_DECREF(module);

		// set the feature names
		std::vector<std::string> featureNames = { "ZScore", "gcContent", "size", "observedCount",
				"biasCorrectedObservedCount", "expectedCountPusm" };
		// convert the feature name vector into python objects
		PyObject* pyFeatureNames = PyList_New(featureNames.size());
		for (size_t i = 0; i < featureNames.size(); ++i) {
			PyList_SetItem(pyFeatureNames, i, PyString_FromString(featureNames[i].c_str()));
		}
		PyObject* pyResultCurrent = PyObject_CallMethod(mlClassifier, "set_features", "O", pyFeatureNames);
		if (!pyResultCurrent) {
			throw std::runtime_error("PYTHON: set_features failed");
		}

		// set the class id's
		PyObject* pyClassIds = PyList_New(3);
		for (size_t i = 0; i < 3; ++i) {
			PyList_SetItem(pyClassIds, i, PyInt_FromLong(i));
		}
		pyResultCurrent = PyObject_CallMethod(mlClassifier, "set_classes", "O", pyClassIds);
		if (!pyResultCurrent) {
			throw std::runtime_error("PYTHON: set_classes failed");
		}

	}
}

void KmerClassificationUnit::storeClassifier(const std::string &filename) {
	if (classificationType == KmerClassificationType::CLASSIFICATION_MACHINE_LEARNING) {
		PyObject* res = PyObject_CallMethod(mlClassifier, "store_classifier", "s", filename.c_str());
		if (!res) {
			throw std::runtime_error("PYTHON: store_classifier failed");
		}
	}
}

void KmerClassificationUnit::loadClassifier(const std::string &filename) {
	if (classificationType == KmerClassificationType::CLASSIFICATION_MACHINE_LEARNING) {
		PyObject* pyResultCurrent = PyObject_CallMethod(mlClassifier, "load_classifier", "s", filename.c_str());
		if (!pyResultCurrent) {
			throw std::runtime_error("PYTHON: load_classifier failed");
		}
	}
}

void KmerClassificationUnit::trainClassifier(Dataset &ds, KmerCounter &genomeCounter) {
	if (classificationType == KmerClassificationType::CLASSIFICATION_MACHINE_LEARNING) {
		featureNames = "ZScore;gcContent;size;observedCount;biasCorrectedObservedCount;expectedCountPusm;type";
		outputPath = ds.plotPath + "kmerTrainData.csv";

		std::ifstream test(outputPath);
		if (test.good()) {
			std::cout << "Alread detected k-mer training data. Skipping extraction...\n";
		} else {
			std::ofstream outfile(outputPath);
			if (!outfile.good()) {
				throw std::runtime_error("Could not create file: " + outputPath);
			}
			outfile << featureNames << "\n";

			std::cout << "Extracting training data for k-mer classification...\n";
			if (length(ds.genome) > 0) { // use reference genome for grabbing k-mers
				extractTrainingDataFromReference(ds.genome, biasUnit.getMinKmerSize(), genomeCounter, outfile);
			} else { // use read dataset for grabbing k-mers? Not sure how, though...
				throw std::runtime_error(
						"Cannot train k-mer classification unit, the reference genome is not provided!");
			}
			std::cout << "Finished extracting training data for k-mer classification.\n";
			outfile.close();
		}

		std::cout << "Choosing the best classifier...\n";
		PyObject* pyResultCurrent = PyObject_CallMethod(mlClassifier, "set_csv_file", "s", outputPath.c_str());
		if (!pyResultCurrent) {
			throw std::runtime_error("PYTHON: set_csv_file failed");
		}
		std::cout << "Finished choosing the best classifier.\n";
	}
}

void KmerClassificationUnit::extractTrainingDataFromReference(const seqan::Dna5String &referenceGenome, size_t k,
		KmerCounter &referenceCounter, std::ofstream &outfile) {
	std::unordered_set<std::string> visitedKmers;

	size_t min_progress = 0;
	size_t num_entries = std::min(500000.0, std::pow(4.0, k) / 2);
	size_t oldK = k;
	size_t numKmerResets = 0;
	size_t maxNumKmerResets = 100;

	size_t numFailed = 0;
	size_t maxFailed = 100000;

	// ensure approximately 1/3 trusted or repetitive k-mers; by using the genome.
	while (visitedKmers.size() < num_entries / 3 && numKmerResets < maxNumKmerResets) {
		size_t randomPosition = rand() % (length(referenceGenome) - k);
		std::string randomKmer = "";
		for (size_t i = 0; i < k; ++i) {
			randomKmer += referenceGenome[randomPosition + i];
			if (visitedKmers.find(randomKmer) == visitedKmers.end()) {
				int numGC = 0;
				for (size_t j = 0; j < k; ++j) {
					if (randomKmer[j] == 'G' || randomKmer[j] == 'C') {
						numGC++;
					}
				}
				double gc = (double) numGC / k;
				size_t covGenome = referenceCounter.countKmerNoRC(randomKmer);
				double covObserved = counter.countKmer(randomKmer);
				double covExpected = pusm.expectedCount(randomKmer).first;
				double covBiasCorrected = 1.0 / biasUnit.getBias(randomKmer) * covObserved;
				double zScore = kmerZScore(randomKmer);
				if (covGenome == 0) {
					writeTrainingString(zScore, gc, k, covObserved, covBiasCorrected, covExpected, KmerType::UNTRUSTED,
							outfile);
				} else if (covGenome == 1) {
					writeTrainingString(zScore, gc, k, covObserved, covBiasCorrected, covExpected, KmerType::TRUSTED,
							outfile);
				} else {
					writeTrainingString(zScore, gc, k, covObserved, covBiasCorrected, covExpected, KmerType::REPEAT,
							outfile);
				}
				visitedKmers.insert(randomKmer);
			} else {
				numFailed++;
				if (numFailed >= maxFailed) {
					k += 2;
					if (k > length(referenceGenome) / 2) {
						k = oldK;
						numKmerResets++;
					}
					break;
				}
			}
		}

		double progress = 100.0 * (double) visitedKmers.size() / num_entries;
		if (progress >= min_progress) {
			std::cout << progress << "%\n";
			min_progress++;
		}
	}

	numFailed = 0;
	numKmerResets = 0;

	k = oldK;
	while (visitedKmers.size() < num_entries && numKmerResets < maxNumKmerResets) {
		std::string randomKmer = "";
		for (size_t j = 0; j < k; ++j) {
			int randomValue = rand() % 4;
			switch (randomValue) {
			case 0:
				randomKmer += "A";
				break;
			case 1:
				randomKmer += "C";
				break;
			case 2:
				randomKmer += "G";
				break;
			default:
				randomKmer += "T";
				break;
			}
		}
		if (visitedKmers.find(randomKmer) == visitedKmers.end()) {
			double covObserved = counter.countKmer(randomKmer);
			if (covObserved == 0) {
				numFailed++;
				if (numFailed >= maxFailed) {
					k += 2;
					numKmerResets++;
					if (k > length(referenceGenome) / 2) {
						k = oldK;
					}
				}
				continue;
			}
			numFailed = 0;
			int numGC = 0;
			for (size_t j = 0; j < k; ++j) {
				if (randomKmer[j] == 'G' || randomKmer[j] == 'C') {
					numGC++;
				}
			}
			double gc = (double) numGC / k;
			size_t covGenome = referenceCounter.countKmerNoRC(randomKmer);
			double covExpected = pusm.expectedCount(randomKmer).first;
			double covBiasCorrected = 1.0 / biasUnit.getBias(randomKmer) * covObserved;
			double zScore = kmerZScore(randomKmer);
			if (covGenome == 0) {
				writeTrainingString(zScore, gc, k, covObserved, covBiasCorrected, covExpected, KmerType::UNTRUSTED,
						outfile);
			} else if (covGenome == 1) {
				writeTrainingString(zScore, gc, k, covObserved, covBiasCorrected, covExpected, KmerType::TRUSTED,
						outfile);
			} else {
				writeTrainingString(zScore, gc, k, covObserved, covBiasCorrected, covExpected, KmerType::REPEAT,
						outfile);
			}
			visitedKmers.insert(randomKmer);
		}

		double progress = 100.0 * (double) visitedKmers.size() / num_entries;
		if (progress >= min_progress) {
			std::cout << progress << "%\n";
			min_progress++;
		}
	}
}

void KmerClassificationUnit::writeTrainingString(double zScore, double gc, size_t k, double countObserved,
		double countBiasCorrected, double countExpectedPusm, KmerType type, std::ofstream &outfile) {
	outfile << zScore << ";" << gc << ";" << k << ";" << countObserved << ";" << countBiasCorrected << ";"
			<< countExpectedPusm << ";" << kmerTypeToNumber(type) << "\n";
}

KmerType KmerClassificationUnit::classifyZScore(double zScore) {
	if (zScore < -2) {
		return KmerType::UNTRUSTED;
	} else if (zScore > 2) {
		return KmerType::REPEAT;
	} else {
		return KmerType::TRUSTED;
	}
}

double gcContent(const std::string &kmer) {
	int num = 0;
	for (size_t i = 0; i < kmer.size(); ++i) {
		if (kmer[i] == 'G' || kmer[i] == 'C') {
			num++;
		}
	}
	return (double) num / kmer.size();
}

KmerType KmerClassificationUnit::classifyKmer(const std::string &kmer) {
	// check if the k-mer is invalid
	if (kmer.find("_") != std::string::npos) {
		throw std::runtime_error("K-mer classification called with an invalid k-mer!");
	}

	if (cachedClassifications.find(kmer) != cachedClassifications.end()) {
		return cachedClassifications[kmer];
	}

	if (classificationType == KmerClassificationType::CLASSIFICATION_STATISTICAL) {
		KmerType type = classifyZScore(kmerZScore(kmer));
		cachedClassifications[kmer] = type;
		return type;
	} else if (classificationType == KmerClassificationType::CLASSIFICATION_NAIVE) {
		double bias = biasUnit.getBias(kmer);
		size_t observedCount = counter.countKmer(kmer);
		double observedCountBiasCorrected = (1 / bias) * observedCount;
		std::pair<double, double> expected = pusm.expectedCount(kmer);
		double expectedCountUnique = expected.first;
		KmerType type;
		if (observedCountBiasCorrected < 0.5 * expectedCountUnique) {
			type = KmerType::UNTRUSTED;
		} else if (observedCountBiasCorrected < 1.5 * expectedCountUnique) {
			type = KmerType::TRUSTED;
		} else {
			type = KmerType::REPEAT;
		}
		cachedClassifications[kmer] = type;
		return type;
	} else if (classificationType == KmerClassificationType::CLASSIFICATION_MACHINE_LEARNING) {
		// build feature vector
		std::vector<double> features;
		features.push_back(kmerZScore(kmer));
		features.push_back(gcContent(kmer));
		features.push_back(kmer.size());
		size_t observedCount = counter.countKmer(kmer);
		features.push_back(observedCount);
		double bias = biasUnit.getBias(kmer);
		double observedCountBiasCorrected = (1 / bias) * observedCount;
		features.push_back(observedCountBiasCorrected);
		std::pair<double, double> expected = pusm.expectedCount(kmer);
		double expectedCountUnique = expected.first;
		features.push_back(expectedCountUnique);
		// convert the vector into Python object
		PyObject* pyFeatures = PyList_New(features.size());
		for (size_t i = 0; i < features.size(); ++i) {
			PyList_SetItem(pyFeatures, i, PyFloat_FromDouble(features[i]));
		}

		/*PyGILState_STATE gstate;
		 gstate = PyGILState_Ensure();*/

		int typeAsInt = -1;

		PyObject* pyResult = PyObject_CallMethod(mlClassifier, "classify", "O", pyFeatures);
		if (!pyResult) {
			throw std::runtime_error("PYTHON: classify failed. kmer was: " + kmer);
		}
		typeAsInt = PyInt_AsLong(pyResult);

		//PyGILState_Release(gstate);

		KmerType kmerType = kmerTypeFromNumber(typeAsInt);
		cachedClassifications[kmer] = kmerType;

		return kmerType;
	} else {
		throw std::runtime_error("Unknown classification type");
	}
}

double KmerClassificationUnit::kmerZScore(const std::string &kmer) {
	// check if the k-mer is invalid
	if (kmer.find("_") != std::string::npos) {
		throw std::runtime_error("K-mer Z-score called with an invalid k-mer!");
	}
	double bias = biasUnit.getBias(kmer);
	size_t observedCount = counter.countKmer(kmer);
	double observedCountBiasCorrected = (1 / bias) * observedCount;
	std::pair<double, double> expected = pusm.expectedCount(kmer);
	double z = (((double) observedCountBiasCorrected) - expected.first) / expected.second;
	if (z != z) {
		throw std::runtime_error(
				"the z-score is not a number!!! observedCount: " + std::to_string(observedCount) + ", bias: "
						+ std::to_string(bias) + ", observedCountBiasCorrected: "
						+ std::to_string(observedCountBiasCorrected) + ", expected.first: "
						+ std::to_string(expected.first) + ", expected.second: " + std::to_string(expected.second));
	}
	return z;
}

size_t KmerClassificationUnit::getMinKmerSize() {
	return biasUnit.getMinKmerSize();
}

KmerClassificationUnit::~KmerClassificationUnit() {
	if (mlClassifier != NULL) {
		Py_DECREF(mlClassifier);
	}
}

void KmerClassificationUnit::clearCache() {
	cachedClassifications.clear();
}
