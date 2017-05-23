/*
 * OverallErrorProfile.h
 *
 *  Created on: Jan 31, 2017
 *      Author: sarah
 */

#pragma once

#include <stddef.h>
#include <map>
#include <string>
#include <unordered_map>
#include <utility>
#include "../external/cereal/types/unordered_map.hpp"
#include "../external/cereal/types/utility.hpp"

#include "../CorrectedRead.h"
#include "../ErrorType.h"
#include "../UtitilyFunctions.hpp"
#include "ErrorProfileUnit.hpp"

class OverallErrorProfile : public ErrorProfileUnit {
public:
	OverallErrorProfile();
	virtual std::unordered_map<ErrorType, double> getErrorProbabilities(const FASTQRead &read, size_t positionInRead);
	virtual std::unordered_map<ErrorType, double> getKmerErrorProbabilities(const std::string &kmer, size_t positionInKmer);
	virtual void loadErrorProfile(const std::string &filepath, KmerCounter &counter);
	virtual void storeErrorProfile(const std::string &filepath);
	virtual void plotErrorProfile();

	virtual void reset();
	virtual void check(const CorrectedRead &corrRead, double acceptProb = 1.0);
	virtual void checkAligned(const CorrectedReadAligned &corrRead, double acceptProb = 1.0);

	virtual void finalize();

	double getOverallErrorRateCurrentBase();
	double getOverallErrorRateNextGap();

	template<class Archive>
	void serialize(Archive & archive) {
		archive(counts, substitutionMatrix, totalCount, noncorrectBases, deletedBases); // serialize things by passing them to the archive
	}
protected:
	virtual std::unordered_map<ErrorType, double> getErrorProbabilitiesFinalized(const FASTQRead &read, size_t positionInRead);
	virtual std::unordered_map<ErrorType, double> getErrorProbabilitiesFinalized(const std::string &kmer, size_t positionInKmer);
private:
	std::unordered_map<ErrorType, size_t> counts;
	std::unordered_map<ErrorType, double> counts_finalized;
	std::unordered_map<std::pair<char, char>, size_t, pairhash> substitutionMatrix;
	std::unordered_map<std::pair<char, char>, double, pairhash> substitutionMatrix_finalized;
	size_t totalCount;
	size_t noncorrectBases;
	size_t deletedBases;

	double overallErrorProbCurrent;
	double overallErrorProbNext;

	bool finalized;
};
