/*
 * PUSM.cpp
 *
 *  Created on: Mar 15, 2017
 *      Author: sarah
 */

#include "PUSM.h"
#include <iostream>

PerfectUniformSequencingModel::PerfectUniformSequencingModel(
		std::shared_ptr<std::unordered_map<size_t, size_t> > &readLengthsPtr, size_t estimatedGenomeSize,
		GenomeType &type) {
	readLengths = readLengthsPtr;
	genomeSize = estimatedGenomeSize;
	genomeType = type;
}

std::pair<double, double> PerfectUniformSequencingModel::expectedCount(const std::string &kmer) {
	if (genomeType == GenomeType::CIRCULAR) {
		return expectedCountCircular(kmer.size());
	} else {
		return expectedCountLinear(kmer.size());
	}
}

std::pair<double, double> PerfectUniformSequencingModel::expectedCountCircular(size_t k) {
	double expected = 0.0;
	double variance = 0.0;
	if (expectedCircular.find(k) != expectedCircular.end()) {
		return expectedCircular[k];
	}
	for (auto pair : (*readLengths.get())) {
		double l = pair.first;
		double n = pair.second;
		if (l < k)
			continue;
		double p = (double) (l - k + 1) / genomeSize;
		if (p > 1.0 || p != p) {
			throw std::runtime_error("Something is wrong with the probability");
		};
		expected += n * p;
		variance += n * p * (1 - p);
	}
	std::pair<double, double> res = std::make_pair(expected, sqrt(variance));
	expectedCircular[k] = res;
	return res;
}

std::pair<double, double> PerfectUniformSequencingModel::expectedCountLinear(size_t k) {
	double expected = 0.0;
	double variance = 0.0;
	if (expectedLinear.find(k) != expectedLinear.end()) {
		return expectedLinear[k];
	}
	for (auto pair : (*readLengths.get())) {
		double l = pair.first;

		if (l > genomeSize) {
			throw std::runtime_error("read length l = " + std::to_string(l) + " > genome size = " + std::to_string(genomeSize));
		}

		double n = pair.second;
		if (l < k)
			continue;
		double p = (1.0 / ((genomeSize - k + 1) * (genomeSize - l + 1))) * (l - k) * (genomeSize - l);
		if (p > 1.0 || p != p) {
			throw std::runtime_error("Something is wrong with the probability. p = " + std::to_string(p) + ", genomeSize = " + std::to_string(genomeSize) + ", k = " + std::to_string(k) + ", l = " + std::to_string(l));
		};
		expected += n * p;
		variance += n * p * (1 - p);
	}
	std::pair<double, double> res = std::make_pair(expected, sqrt(variance));
	expectedLinear[k] = res;
	return res;
}
