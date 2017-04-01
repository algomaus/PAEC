/*
 * MotifErrorProfile.h
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
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>

#include "../../ErrorType.h"
#include "../../UtitilyFunctions.hpp"
#include "../ErrorProfileUnit.hpp"
#include "../../KmerClassification/KmerCounter.h"
#include "../../CorrectedRead.h"

#include "MotifTree.h"

static const int MAX_MOTIF_SIZE = 6;


class MotifErrorProfile : public ErrorProfileUnit {
public:
	MotifErrorProfile(KmerCounter &kmerCounter);
	virtual std::unordered_map<ErrorType, double> getErrorProbabilities(const FASTQRead &read, size_t positionInRead);
	virtual std::unordered_map<ErrorType, double> getKmerErrorProbabilities(const std::string &kmer, size_t positionInKmer);
	virtual void loadErrorProfile(const std::string &filepath, KmerCounter &counter);
	virtual void storeErrorProfile(const std::string &filepath);
	virtual void plotErrorProfile();

	virtual void reset();
	virtual void check(const CorrectedRead &corrRead, double acceptProb = 1.0);
	virtual void checkAligned(const CorrectedReadAligned &corrRead, double acceptProb = 1.0);

	virtual void finalize();

	double findMostSignificantZScore(const ErrorType &type, const std::string &sequence, int posInSequence);

	template<class Archive>
	void serialize(Archive & archive) {
		archive(motifTree, finalized); // serialize things by passing them to the archive
	}
protected:
	virtual std::unordered_map<ErrorType, double> getErrorProbabilitiesFinalized(const FASTQRead &read, size_t positionInRead);
	virtual std::unordered_map<ErrorType, double> getErrorProbabilitiesFinalized(const std::string &kmer, size_t positionInKmer);
private:
	void updateMotifData(int positionInSequence, const std::string &sequence, const ErrorType &type);
	void computeZScores();

	MotifTree motifTree;
	bool finalized;
	KmerCounter &counter;
};
