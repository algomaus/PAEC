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
	endPos = 0;
}

std::string CorrectedReadAligned::toString() const {
	std::stringstream ss;
	if (correctedRead.sequence.size() > 0) {
		ss << originalRead << "\n" << correctedRead << "\n";
		ss << beginPos; // << " " << positionOffset << " " << alignedCorrections.size();
	}
	return ss.str();
}

CorrectedReadAligned::CorrectedReadAligned(const FASTQRead &original, size_t mappingPos) :
		CorrectedRead(original) {
	beginPos = mappingPos;
	endPos = 0;
}

// TODO: Correct chimeric breaks, too! This is a bit uglier because it means the read has to be splitted.

// Quality scores: Simply leave them as-they-are for substitutions, remove them for insertions and add the average value for deletions
void CorrectedReadAligned::applyCorrection(const CorrectionAligned &ca) {
	Correction corr = ca.correction;
	assert(corr.positionInRead < correctedRead.sequence.size());
	correctedRead.sequence.replace(corr.positionInRead, corr.originalBases.size(), corr.correctedBases);
	alignedCorrections.push_back(ca);

	if (corr.type == ErrorType::INSERTION) {
		correctedRead.quality.replace(corr.positionInRead, 1, "");
		originalPositions.erase(originalPositions.begin() + corr.positionInRead);
	} else if (corr.type == ErrorType::DEL_OF_A || corr.type == ErrorType::DEL_OF_C || corr.type == ErrorType::DEL_OF_G
			|| corr.type == ErrorType::DEL_OF_T) {
		// find average quality score
		char qualLeft = correctedRead.quality[corr.positionInRead];
		char qualRight = correctedRead.quality[corr.positionInRead + 1];
		char qualMiddle = (qualLeft + qualRight) / 2;
		std::string newQual = "";
		newQual += qualLeft;
		newQual += qualMiddle;
		correctedRead.quality.replace(corr.positionInRead, 1, newQual);
		// update offset
		originalPositions.insert(originalPositions.begin() + corr.positionInRead, originalPositions[corr.positionInRead]);
	} else if (corr.type == ErrorType::MULTIDEL && corr.correctedBases.size() > 0) {
		// find average quality score
		char qualLeft = correctedRead.quality[corr.positionInRead];
		char qualRight = correctedRead.quality[corr.positionInRead + 1];
		char qualMiddle = (qualLeft + qualRight) / 2;
		std::string newQual = "";
		newQual += qualLeft;
		for (size_t i = 0; i < corr.correctedBases.size() - 1; ++i) {
			newQual += qualMiddle;
		}
		correctedRead.quality.replace(corr.positionInRead, 1, newQual);
		// update offset
		for (size_t i = 0; i < corr.correctedBases.size(); ++i) {
			originalPositions.insert(originalPositions.begin() + corr.positionInRead,
					originalPositions[corr.positionInRead]);
		}
	}

	assert(correctedRead.sequence.size() == correctedRead.quality.size());
}
