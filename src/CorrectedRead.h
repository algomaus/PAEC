/*
 * CorrectedRead.h
 *
 *  Created on: Jan 31, 2017
 *      Author: sarah
 */

#pragma once

#include <stddef.h>
#include <string>
#include <vector>

#include <cereal/types/vector.hpp>

#include "Correction.h"
#include "ErrorType.h"
#include "FASTQRead.h"

class CorrectedRead {
public:
	CorrectedRead();
	CorrectedRead(const FASTQRead &original);
	void applyCorrection(const Correction &corr);
	void applyCorrection(const ErrorType &type, size_t posInCorrectedRead, double prob = 1);
	FASTQRead getOriginalRead();
	std::string toString() const;

	FASTQRead originalRead;
	FASTQRead correctedRead;
	int positionOffset;
	std::vector<Correction> corrections;

	template<class Archive>
	void serialize(Archive & archive) {
		archive(originalRead, correctedRead, positionOffset, corrections); // serialize things by passing them to the archive
	}
};

inline std::ostream& operator<<(std::ostream &os, const CorrectedRead &corr) {
	return os << corr.toString();
}

inline std::istream& operator>>(std::istream &is, CorrectedRead &corr) {
	size_t corrSize;
	is >> corr.originalRead >> corr.correctedRead >> corr.positionOffset >> corrSize;
	for (size_t i = 0; i < corrSize; ++i) {
		Correction co;
		is >> co;
		corr.corrections.push_back(co);
	}
	return is;
}
