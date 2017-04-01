/*
 * CorrectedReadAligned.h
 *
 *  Created on: Jan 31, 2017
 *      Author: sarah
 */

#pragma once

#include <memory>

#include <cereal/types/vector.hpp>

#include "../CorrectedRead.h"
#include "CorrectionAligned.h"

class CorrectedReadAligned: public CorrectedRead {
public:
	CorrectedReadAligned();
	CorrectedReadAligned(const FASTQRead &original, size_t mappingPos);
	void applyCorrection(const CorrectionAligned &ca);
	std::string toString() const;
	size_t beginPos;
	std::vector<CorrectionAligned> alignedCorrections;

	template<class Archive>
	void serialize(Archive & archive) {
		archive(originalRead, correctedRead, positionOffset, beginPos, alignedCorrections); // serialize things by passing them to the archive
	}
};

inline std::ostream& operator<<(std::ostream &os, const CorrectedReadAligned &corr) {
	return os << corr.toString();
}

inline std::istream& operator>>(std::istream &is, CorrectedReadAligned &corr) {
	size_t corrSize;
	is >> corr.originalRead >> corr.correctedRead >> corr.beginPos >> corr.positionOffset >> corrSize;
	for (size_t i = 0; i < corrSize; ++i) {
		CorrectionAligned ca;
		is >> ca;
		corr.alignedCorrections.push_back(ca);
	}
	return is;
}

