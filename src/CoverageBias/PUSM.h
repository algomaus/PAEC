/*
 * PUSM.h
 *
 *  Created on: Feb 6, 2017
 *      Author: sarah
 */

#pragma once

#include <stddef.h>
#include <cmath>
#include <string>
#include <unordered_map>
#include <memory>
#include <utility>

#include "GenomeType.h"

class PerfectUniformSequencingModel {
public:
	PerfectUniformSequencingModel(std::shared_ptr<std::unordered_map<size_t, size_t> > &readLengthsPtr,
			size_t estimatedGenomeSize, GenomeType &type);
	std::pair<double, double> expectedCount(const std::string &kmer);
private:
	std::pair<double, double> expectedCountLinear(size_t k);
	std::pair<double, double> expectedCountCircular(size_t k);

	std::shared_ptr<std::unordered_map<size_t, size_t> > readLengths;
	size_t genomeSize;
	GenomeType genomeType;
	std::unordered_map<size_t, std::pair<double, double> > expectedLinear;
	std::unordered_map<size_t, std::pair<double, double> > expectedCircular;
};

