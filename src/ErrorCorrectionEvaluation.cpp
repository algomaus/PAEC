
#include "ErrorCorrectionEvaluation.h"

ErrorCorrectionEvaluation::ErrorCorrectionEvaluation() {
	totalTruePositives = 0;
	totalFalsePositives = 0;
	totalTrueNegatives = 0;
	totalFalseNegatives = 0;
	for (ErrorType type : errorTypeIterator()) {
		truePositives[type] = 0;
		falsePositives[type] = 0;
		trueNegatives[type] = 0;
		falseNegatives[type] = 0;
	}
}

void ErrorCorrectionEvaluation::checkAligned(const CorrectedReadAligned &cra) {
	trueReads.push(cra);
}

void ErrorCorrectionEvaluation::finalize() {
	std::cout << "totalTruePositives: " << totalTruePositives << "\n";
	std::cout << "totalFalsePositives: " << totalFalsePositives << "\n";
	std::cout << "totalTrueNegatives: " << totalTrueNegatives << "\n";
	std::cout << "totalFalseNegatives: " << totalFalseNegatives << "\n";

	double totalGain = (totalTruePositives - totalFalsePositives) / (double) (totalTruePositives + totalFalseNegatives);
	double totalSensitivity = (totalTruePositives) / (double) (totalTruePositives + totalFalseNegatives);
	double totalSpecificity = (totalTrueNegatives) / (double) (totalTrueNegatives + totalFalsePositives);

	std::cout << "totalGain: " << totalGain << "\n";
	std::cout << "totalSensitivity: " << totalSensitivity << "\n";
	std::cout << "totalSpecificity: " << totalSpecificity << "\n";
	double fscoreTotal = (2 * totalTruePositives)
			/ (double) (2 * totalTruePositives + totalFalsePositives + totalFalseNegatives);
	std::cout << "fscoreTotal: " << fscoreTotal << "\n";

	size_t truePositivesSubstitutions = truePositives[ErrorType::SUB_FROM_A] + truePositives[ErrorType::SUB_FROM_C]
			+ truePositives[ErrorType::SUB_FROM_G] + truePositives[ErrorType::SUB_FROM_T];
	size_t falsePositivesSubstitutions = falsePositives[ErrorType::SUB_FROM_A] + falsePositives[ErrorType::SUB_FROM_C]
			+ falsePositives[ErrorType::SUB_FROM_G] + falsePositives[ErrorType::SUB_FROM_T];
	size_t trueNegativesSubstitutions = trueNegatives[ErrorType::SUB_FROM_A] + trueNegatives[ErrorType::SUB_FROM_C]
			+ trueNegatives[ErrorType::SUB_FROM_G] + trueNegatives[ErrorType::SUB_FROM_T];
	size_t falseNegativesSubstitutions = falseNegatives[ErrorType::SUB_FROM_A] + falseNegatives[ErrorType::SUB_FROM_C]
			+ falseNegatives[ErrorType::SUB_FROM_G] + falseNegatives[ErrorType::SUB_FROM_T];
	double substutitionsGain = (truePositivesSubstitutions - falsePositivesSubstitutions)
			/ (double) (truePositivesSubstitutions + falseNegativesSubstitutions);
	double substutitionsSensitivity = (truePositivesSubstitutions)
			/ (double) (truePositivesSubstitutions + falseNegativesSubstitutions);
	double substutitionsSpecificity = (trueNegativesSubstitutions)
			/ (double) (trueNegativesSubstitutions + falsePositivesSubstitutions);
	std::cout << "truePositivesSubstitutions: " << truePositivesSubstitutions << "\n";
	std::cout << "falsePositivesSubstitutions: " << falsePositivesSubstitutions << "\n";
	std::cout << "trueNegativesSubstitutions: " << trueNegativesSubstitutions << "\n";
	std::cout << "falseNegativesSubstitutions: " << falseNegativesSubstitutions << "\n";
	std::cout << "substutitionsGain: " << substutitionsGain << "\n";
	std::cout << "substutitionsSensitivity: " << substutitionsSensitivity << "\n";
	std::cout << "substutitionsSpecificity: " << substutitionsSpecificity << "\n";
	double fscoreSubstitutions = (2 * truePositivesSubstitutions)
			/ (double) (2 * truePositivesSubstitutions + falsePositivesSubstitutions + falseNegativesSubstitutions);
	std::cout << "fscoreSubstitutions: " << fscoreSubstitutions << "\n";

	size_t truePositivesInsertions = truePositives[ErrorType::INSERTION];
	size_t falsePositivesInsertions = falsePositives[ErrorType::INSERTION];
	size_t trueNegativesInsertions = trueNegatives[ErrorType::INSERTION];
	size_t falseNegativesInsertions = falseNegatives[ErrorType::INSERTION];
	double insertionsGain = (truePositivesInsertions - falsePositivesInsertions)
			/ (double) (truePositivesInsertions + falseNegativesInsertions);
	double insertionsSensitivity = (truePositivesInsertions)
			/ (double) (truePositivesInsertions + falseNegativesInsertions);
	double insertionsSpecificity = (trueNegativesInsertions)
			/ (double) (trueNegativesInsertions + falsePositivesInsertions);
	std::cout << "truePositivesInsertions: " << truePositivesInsertions << "\n";
	std::cout << "falsePositivesInsertions: " << falsePositivesInsertions << "\n";
	std::cout << "trueNegativesInsertions: " << trueNegativesInsertions << "\n";
	std::cout << "falseNegativesInsertions: " << falseNegativesInsertions << "\n";
	std::cout << "insertionsGain: " << insertionsGain << "\n";
	std::cout << "insertionsSensitivity: " << insertionsSensitivity << "\n";
	std::cout << "insertionsSpecificity: " << insertionsSpecificity << "\n";
	double fscoreInsertions = (2 * truePositivesInsertions)
			/ (double) (2 * truePositivesInsertions + falsePositivesInsertions + falseNegativesInsertions);
	std::cout << "fscoreInsertions: " << fscoreInsertions << "\n";

	size_t truePositivesDeletions = truePositives[ErrorType::DEL_OF_A] + truePositives[ErrorType::DEL_OF_C]
			+ truePositives[ErrorType::DEL_OF_G] + truePositives[ErrorType::DEL_OF_T]
			+ truePositives[ErrorType::MULTIDEL];
	size_t falsePositivesDeletions = falsePositives[ErrorType::DEL_OF_A] + falsePositives[ErrorType::DEL_OF_C]
			+ falsePositives[ErrorType::DEL_OF_G] + falsePositives[ErrorType::DEL_OF_T]
			+ falsePositives[ErrorType::MULTIDEL];
	size_t trueNegativesDeletions = trueNegatives[ErrorType::DEL_OF_A] + trueNegatives[ErrorType::DEL_OF_C]
			+ trueNegatives[ErrorType::DEL_OF_G] + trueNegatives[ErrorType::DEL_OF_T]
			+ trueNegatives[ErrorType::MULTIDEL];
	size_t falseNegativesDeletions = falseNegatives[ErrorType::DEL_OF_A] + falseNegatives[ErrorType::DEL_OF_C]
			+ falseNegatives[ErrorType::DEL_OF_G] + falseNegatives[ErrorType::DEL_OF_T]
			+ falseNegatives[ErrorType::MULTIDEL];
	double deletionsGain = (truePositivesDeletions - falsePositivesDeletions)
			/ (double) (truePositivesDeletions + falseNegativesDeletions);
	double deletionsSensitivity = (truePositivesDeletions)
			/ (double) (truePositivesDeletions + falseNegativesDeletions);
	double deletionsSpecificity = (trueNegativesDeletions)
			/ (double) (trueNegativesDeletions + falsePositivesDeletions);
	std::cout << "truePositivesDeletions: " << truePositivesDeletions << "\n";
	std::cout << "falsePositivesDeletions: " << falsePositivesDeletions << "\n";
	std::cout << "trueNegativesDeletions: " << trueNegativesDeletions << "\n";
	std::cout << "falseNegativesDeletions: " << falseNegativesDeletions << "\n";
	std::cout << "deletionsGain: " << deletionsGain << "\n";
	std::cout << "deletionsSensitivity: " << deletionsSensitivity << "\n";
	std::cout << "deletionsSpecificity: " << deletionsSpecificity << "\n";
	double fscoreDeletions = (2 * truePositivesDeletions)
			/ (double) (2 * truePositivesDeletions + falsePositivesDeletions + falseNegativesDeletions);
	std::cout << "fscoreDeletions: " << fscoreDeletions << "\n";

}

void ErrorCorrectionEvaluation::check(const CorrectedRead &cr) {
	// find the corresponding true read
	if (trueReads.front().correctedRead.id < cr.correctedRead.id) {
		//break; // since the corrected reads arrive in a sorted order
	} else if (trueReads.front().correctedRead.id == cr.correctedRead.id) {
		// found the corresponding true read
		fillConfusionMatrix(trueReads.front(), cr);
		trueReads.pop();
	}
}

void ErrorCorrectionEvaluation::fillConfusionMatrix(const CorrectedReadAligned &truth, const CorrectedRead &cr) {
	assert(truth.originalRead.sequence.size() == cr.originalRead.sequence.size());

	std::vector<bool> unalteredBasesTruth(truth.originalRead.sequence.size(), true);
	std::vector<bool> unalteredBasesRead(cr.originalRead.sequence.size(), true);

	for (CorrectionAligned corrAl: truth.alignedCorrections) {
		unalteredBasesTruth[corrAl.correction.originalReadPos] = false;
		bool found = false;
		for (Correction corr : cr.corrections) {
			if (corrAl.correction.originalReadPos == corr.originalReadPos && corrAl.correction.type == corr.type) {
				found = true;
				break;
			}
		}
		if (found) {
			truePositives[corrAl.correction.type]++;
			totalTruePositives++;
		} else {
			falseNegatives[corrAl.correction.type]++;
			totalFalseNegatives++;
		}
	}

	for (Correction corr : cr.corrections) {
		unalteredBasesRead[corr.originalReadPos] = false;
		bool found = false;
		for (CorrectionAligned corrAl: truth.alignedCorrections) {
			if (corrAl.correction.originalReadPos == corr.originalReadPos && corrAl.correction.type == corr.type) {
				found = true;
				break;
			}
		}

		if (!found) {
			if (corr.originalReadPos >= truth.correctedRead.sequence.size() || truth.correctedRead.sequence[corr.originalReadPos] != 'S') { // ignore soft-clipped bases
				falsePositives[corr.type]++;
				totalFalsePositives++;
			}
		}
	}

	for (size_t i = 0; i < unalteredBasesTruth.size(); ++i) {
		if (i < truth.correctedRead.sequence.size() && truth.correctedRead.sequence[i] == 'S') {
			continue;
		}
		if (unalteredBasesTruth[i] && unalteredBasesRead[i]) {
			totalTrueNegatives++;
			if (i < unalteredBasesTruth.size() - 1) {
				totalTrueNegatives++; // also consider true negative gaps (potential positions of deletions)
			}
			for (ErrorType type : errorTypeIterator()) {
				if (type == ErrorType::INSERTION) {
					trueNegatives[type]++;
				} else if (type == ErrorType::SUB_FROM_A || type == ErrorType::SUB_FROM_C || type == ErrorType::SUB_FROM_G || type == ErrorType::SUB_FROM_T) {
					trueNegatives[type]++;
				} else if (type == ErrorType::DEL_OF_A || type == ErrorType::DEL_OF_C || type == ErrorType::DEL_OF_G || type == ErrorType::DEL_OF_T) {
					if (i < unalteredBasesTruth.size() - 1) { // since no deletions can happen at the last base of a read
						trueNegatives[type]++;
					}
				}
			}
		}
	}
}
