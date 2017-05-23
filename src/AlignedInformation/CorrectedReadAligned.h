/*
 * CorrectedReadAligned.h
 *
 *  Created on: Jan 31, 2017
 *      Author: sarah
 */

#pragma once

#include <memory>

#include "../external/cereal/types/vector.hpp"

#include "../CorrectedRead.h"
#include "CorrectionAligned.h"

class CorrectedReadAligned: public CorrectedRead {
public:
	CorrectedReadAligned();
	CorrectedReadAligned(const FASTQRead &original, size_t mappingPos);
	void applyCorrection(const CorrectionAligned &ca);
	std::string toString() const;
	size_t beginPos;
	size_t endPos;
	std::vector<CorrectionAligned> alignedCorrections;

	template<class Archive>
	void serialize(Archive & archive) {
		archive(originalRead, correctedRead, originalPositions, beginPos, endPos, alignedCorrections); // serialize things by passing them to the archive
	}
};

