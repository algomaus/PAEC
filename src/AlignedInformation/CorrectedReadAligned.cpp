/*
 * CorrectedReadAligned.cpp
 *
 *  Created on: Jan 31, 2017
 *      Author: sarah
 */

#include "CorrectedReadAligned.h"

#include <cassert>
#include <sstream>
#include <string>
#include <vector>

CorrectedReadAligned::CorrectedReadAligned() :
		CorrectedRead() {
	beginPos = 0;
}

std::string CorrectedReadAligned::toString() const {
	std::stringstream ss;
	if (correctedRead.sequence.size() > 0) {
		ss << originalRead << "\n" << correctedRead << "\n";
		ss << beginPos << " " << positionOffset << " " << alignedCorrections.size();
		for (CorrectionAligned ca : alignedCorrections) {
			ss << "\n" << ca;
		}
	}
	return ss.str();
}

CorrectedReadAligned::CorrectedReadAligned(const FASTQRead &original, size_t mappingPos) :
		CorrectedRead(original) {
	beginPos = mappingPos;
}

// TODO: Correct chimeric breaks, too! This is a bit uglier because it means the read has to be splitted.

// Quality scores: Simply leave them as-they-are for substitutions, remove them for insertions and add the average value for deletions
void CorrectedReadAligned::applyCorrection(const CorrectionAligned &ca) {
	Correction corr = ca.correction;
	assert(corr.positionInRead + positionOffset < correctedRead.sequence.size());

	correctedRead.sequence.replace(corr.positionInRead + positionOffset, corr.originalBases.size(), corr.correctedBases);
	alignedCorrections.push_back(ca);

	if (corr.type == ErrorType::INSERTION) {
		correctedRead.quality.replace(corr.positionInRead + positionOffset, 1, "");
		positionOffset--;
	} else if (corr.type == ErrorType::DEL_OF_A || corr.type == ErrorType::DEL_OF_C || corr.type == ErrorType::DEL_OF_G
			|| corr.type == ErrorType::DEL_OF_T) {
		// find average quality score
		char qualLeft = correctedRead.quality[corr.positionInRead + positionOffset];
		char qualRight = correctedRead.quality[corr.positionInRead + positionOffset + 1];
		char qualMiddle = (qualLeft + qualRight) / 2;
		std::string newQual = "";
		newQual += qualLeft;
		newQual += qualMiddle;
		correctedRead.quality.replace(corr.positionInRead + positionOffset, 1, newQual);
		// update offset
		positionOffset++;
	} else if (corr.type == ErrorType::MULTIDEL && corr.correctedBases.size() > 0) {
		// find average quality score
		char qualLeft = correctedRead.quality[corr.positionInRead + positionOffset];
		char qualRight = correctedRead.quality[corr.positionInRead + positionOffset + 1];
		char qualMiddle = (qualLeft + qualRight) / 2;
		std::string newQual = "";
		newQual += qualLeft;
		for (size_t i = 0; i < corr.correctedBases.size() - 1; ++i) {
			newQual += qualMiddle;
		}
		correctedRead.quality.replace(corr.positionInRead + positionOffset, 1, newQual);
		// update offset
		positionOffset += corr.correctedBases.size() - 1;
	}

	assert(correctedRead.sequence.size() == correctedRead.quality.size());
}
