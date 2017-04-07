#pragma once

#include <queue>
#include <cassert>
#include <unordered_map>
#include "ErrorType.h"
#include "AlignedInformation/CorrectedReadAligned.h"
#include "CorrectedRead.h"

class ErrorCorrectionEvaluation {
public:
	ErrorCorrectionEvaluation();
	void check(const CorrectedRead &corrRead);
	void checkAligned(const CorrectedReadAligned &corrRead);
	void finalize();
private:
	void fillConfusionMatrix(const CorrectedReadAligned &truth, const CorrectedRead &cr);

	std::queue<CorrectedReadAligned> trueReads;

	std::unordered_map<ErrorType, size_t> truePositives;
	std::unordered_map<ErrorType, size_t> falsePositives;
	std::unordered_map<ErrorType, size_t> trueNegatives;
	std::unordered_map<ErrorType, size_t> falseNegatives;
	size_t totalTruePositives;
	size_t totalFalsePositives;
	size_t totalTrueNegatives;
	size_t totalFalseNegatives;
};
