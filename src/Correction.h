/*
 * Correction.h
 *
 *  Created on: Jan 31, 2017
 *      Author: sarah
 */

#pragma once

#include <stddef.h>
#include <iostream>
#include <string>
#include <cereal/types/string.hpp>

#include "ErrorType.h"

class Correction {
public:
	Correction();
	Correction(size_t pos, size_t origPos, const std::string &from, const std::string &to, double prob);
	Correction(size_t pos, size_t origPos, const std::string &from, const std::string &to);
	Correction(size_t pos, size_t origPos, const std::string &from, const std::string &to, const ErrorType &errorType);
	Correction(size_t pos, size_t origPos, const std::string &from, const std::string &to, double prob, const ErrorType &errorType);
	std::string toString() const;

	size_t positionInRead; // position in the corrected read
	std::string originalBases;
	std::string correctedBases;
	double correctionProbability;
	ErrorType type;
	size_t originalReadPos; // position in the original read

	ErrorType inferErrorType(const std::string &from, const std::string &to);

	template<class Archive>
	void serialize(Archive & archive) {
		archive(positionInRead, originalReadPos, originalBases, correctedBases, correctionProbability, type); // serialize things by passing them to the archive
	}
};

inline std::ostream& operator<<(std::ostream & os, const Correction &co) {
	return os << co.toString();
}

inline std::istream& operator>>(std::istream & is, Correction &co) {
	is >> co.type >> co.positionInRead >> co.originalBases >> co.correctedBases >> co.correctionProbability;

	if (co.originalBases == "_") {
		co.originalBases = "";
	}
	if (co.correctedBases == "_") {
		co.correctedBases = "";
	}

	return is;
}
