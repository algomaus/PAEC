/*
 * CorrectedRead.cpp
 *
 *  Created on: Jan 31, 2017
 *      Author: sarah
 */

#include <cassert>
#include <sstream>
#include "CorrectedRead.h"


CorrectedRead::CorrectedRead() {
	correctedRead = FASTQRead();
	positionOffset = 0;
	originalRead = FASTQRead();
}

std::string CorrectedRead::toString() const {
	std::stringstream ss;
	if (corrections.size() > 0) {
		ss << originalRead << "\n" << correctedRead << "\n";
		ss << positionOffset << " " << corrections.size();
		for (Correction co : corrections) {
			ss << "\n" << co;
		}
	}
	return ss.str();
}

CorrectedRead::CorrectedRead(const FASTQRead &original) {
	correctedRead = original;
	positionOffset = 0;
	originalRead = original;
	assert(correctedRead.sequence.size() == correctedRead.quality.size());
}

FASTQRead CorrectedRead::getOriginalRead() {
	return originalRead;
}

void CorrectedRead::applyCorrection(const ErrorType &type, size_t posInCorrectedRead, double prob) {
	std::string fromBases = "";
	fromBases += correctedRead.sequence[posInCorrectedRead];
	//std::cout << "correctedRead.sequence[posInCorrectedRead]: " << correctedRead.sequence[posInCorrectedRead] << "\n";
	std::string toBases = "";
	if (type == ErrorType::INSERTION) {
		applyCorrection(Correction(posInCorrectedRead - positionOffset, fromBases, toBases, ErrorType::INSERTION));
	} else if (type == ErrorType::SUB_FROM_A) {
		toBases += 'A';
		applyCorrection(Correction(posInCorrectedRead - positionOffset, fromBases, toBases, ErrorType::SUB_FROM_A));
	} else if (type == ErrorType::SUB_FROM_C) {
		toBases += 'C';
		applyCorrection(Correction(posInCorrectedRead - positionOffset, fromBases, toBases, ErrorType::SUB_FROM_C));
	} else if (type == ErrorType::SUB_FROM_G) {
		toBases += 'G';
		applyCorrection(Correction(posInCorrectedRead - positionOffset, fromBases, toBases, ErrorType::SUB_FROM_G));
	} else if (type == ErrorType::SUB_FROM_T) {
		toBases += 'T';
		applyCorrection(Correction(posInCorrectedRead - positionOffset, fromBases, toBases, ErrorType::SUB_FROM_T));
	} else {
		toBases += fromBases;
		if (type == ErrorType::DEL_OF_A) {
			toBases += 'A';
			applyCorrection(Correction(posInCorrectedRead - positionOffset, fromBases, toBases, ErrorType::DEL_OF_A));
		} else if (type == ErrorType::DEL_OF_C) {
			toBases += 'C';
			applyCorrection(Correction(posInCorrectedRead - positionOffset, fromBases, toBases, ErrorType::DEL_OF_C));
		} else if (type == ErrorType::DEL_OF_G) {
			toBases += 'G';
			applyCorrection(Correction(posInCorrectedRead - positionOffset, fromBases, toBases, ErrorType::DEL_OF_G));
		} else if (type == ErrorType::DEL_OF_T) {
			toBases += 'T';
			applyCorrection(Correction(posInCorrectedRead - positionOffset, fromBases, toBases, ErrorType::DEL_OF_T));
		} else { // deletion of multiple bases or a chimeric break
			toBases += '_';
			applyCorrection(Correction(posInCorrectedRead - positionOffset, fromBases, toBases, ErrorType::MULTIDEL));
		}
	}
	assert(correctedRead.sequence.size() == correctedRead.quality.size());
}

// Quality scores: Simply leave them as-they-are for substitutions, remove them for insertions and add the average value for deletions
void CorrectedRead::applyCorrection(const Correction &corr) {
	assert(corr.originalBases != corr.correctedBases);

	//std::cout << "corr.originalBases: " << corr.originalBases << "\n";
	//std::cout << "corr.correctedBases: " << corr.correctedBases << "\n";

	correctedRead.sequence.replace(corr.positionInRead + positionOffset, corr.originalBases.size(), corr.correctedBases);
	corrections.push_back(corr);

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

	//std::cout << "corr.type: " << corr.type << "\n";
	//std::cout << "correctedRead.sequence.size(): " << correctedRead.sequence.size() << "\n";
	//std::cout << "correctedRead.quality.size(): " << correctedRead.quality.size() << "\n";
	assert(correctedRead.sequence.size() == correctedRead.quality.size());
}
