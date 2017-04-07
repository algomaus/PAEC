/*
 * Correction.cpp
 *
 *  Created on: Jan 31, 2017
 *      Author: sarah
 */

#include "Correction.h"

#include <sstream>
#include <stdexcept>
#include <cassert>

Correction::Correction() {
	positionInRead = 0;
	originalBases = "";
	correctedBases = "";
	correctionProbability = 1;
	type = ErrorType::CORRECT;
	originalReadPos = 0;
}

std::string Correction::toString() const {
	std::stringstream ss;
	ss << type << " " << positionInRead << " ";
	if (originalBases.empty()) {
		ss << ".";
	} else {
		ss << originalBases;
	}
	ss << " ";
	if (correctedBases.empty()) {
		ss << ".";
	} else {
		ss << correctedBases;
	}
	ss << " " << correctionProbability;
	return ss.str();
}

Correction::Correction(size_t pos, size_t origPos, const std::string &from, const std::string &to, double prob) {
	assert(from != to);
	positionInRead = pos;
	originalBases = from;
	correctedBases = to;
	correctionProbability = prob;
	type = inferErrorType(from, to);
	originalReadPos = origPos;
}

Correction::Correction(size_t pos, size_t origPos, const std::string &from, const std::string &to, double prob, const ErrorType &errorType) {
	assert(from != to);
	positionInRead = pos;
	originalBases = from;
	correctedBases = to;
	correctionProbability = prob;
	type = errorType;
	originalReadPos = origPos;
}

Correction::Correction(size_t pos, size_t origPos, const std::string &from, const std::string &to, const ErrorType &errorType) {
	assert(errorType == ErrorType::MULTIDEL || from != to);
	positionInRead = pos;
	originalBases = from;
	correctedBases = to;
	correctionProbability = 1;
	type = errorType;
	originalReadPos = origPos;
}

Correction::Correction(size_t pos, size_t origPos, const std::string &from, const std::string &to) {
	assert(from != to);
	positionInRead = pos;
	originalBases = from;
	correctedBases = to;
	correctionProbability = 1;
	type = inferErrorType(from, to);
	originalReadPos = origPos;
}

ErrorType Correction::inferErrorType(const std::string &from, const std::string &to) {
	assert(from != to);
	if (from.size() != 1) {
		throw std::runtime_error("from.size() != 0");
	}

	ErrorType type;

	if (to.size() == 0) {
			type = ErrorType::INSERTION;
	} else if (to.size() == 1) {
		if (to[0] == 'A') {
			type = ErrorType::SUB_FROM_A;
		} else if (to[0] == 'C') {
			type = ErrorType::SUB_FROM_C;
		} else if (to[0] == 'G') {
			type = ErrorType::SUB_FROM_G;
		} else {
			type = ErrorType::SUB_FROM_T;
		}
	} else if (to.size() == 2) { // Deletion of A, C, G, T
		if (to[1] == 'A') {
			type = ErrorType::DEL_OF_A;
		} else if (to[1] == 'C') {
			type = ErrorType::DEL_OF_C;
		} else if (to[1] == 'G') {
			type = ErrorType::DEL_OF_G;
		} else {
			type = ErrorType::DEL_OF_T;
		}
	} else if (to.size() > 2) {
		type = ErrorType::MULTIDEL;
	}

	return type;
}
