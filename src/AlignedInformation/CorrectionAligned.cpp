/*
 * CorrectionAligned.cpp
 *
 *  Created on: Jan 31, 2017
 *      Author: sarah
 */

#include "CorrectionAligned.h"

CorrectionAligned::CorrectionAligned() {
	positionInReference = 0;
	correction = Correction();
}

CorrectionAligned::CorrectionAligned(size_t posInRef, const Correction &corr) {
	positionInReference = posInRef;
	correction = corr;
}
