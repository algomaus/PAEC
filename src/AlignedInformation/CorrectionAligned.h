/*
 * CorrectionAligned.h
 *
 *  Created on: Jan 31, 2017
 *      Author: sarah
 */

#pragma once

#include "../Correction.h"

class CorrectionAligned {
public:
	CorrectionAligned();

	CorrectionAligned(size_t posInRef, const Correction &corr);

	Correction correction;
	size_t positionInReference;

	template<class Archive>
	void serialize(Archive & archive) {
		archive(correction, positionInReference); // serialize things by passing them to the archive
	}
};
