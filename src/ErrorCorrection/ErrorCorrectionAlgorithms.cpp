#include <stddef.h>
#include <cassert>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <algorithm>

#include "../CorrectedRead.h"
#include "../ErrorProfile/ErrorProfileUnit.hpp"
#include "../ErrorType.h"
#include "../FASTQRead.h"
#include "../KmerClassification/KmerClassificationUnit.h"
#include "../KmerClassification/KmerType.h"

std::pair<ErrorType, double> mostLikelyCurrentBase(std::unordered_map<ErrorType, double> &errorProbabilities) {
	double bestProb = errorProbabilities[ErrorType::CORRECT];
	ErrorType bestType = ErrorType::CORRECT;
	for (ErrorType type : errorTypesCurrentBase()) {
		if (errorProbabilities[type] > bestProb) {
			bestProb = errorProbabilities[type];
			bestType = type;
		}
	}
	return std::make_pair(bestType, bestProb);
}

std::pair<ErrorType, double> mostLikelyNextGap(std::unordered_map<ErrorType, double> &errorProbabilities) {
	double bestProb = errorProbabilities[ErrorType::NODEL];
	ErrorType bestType = ErrorType::NODEL;
	for (ErrorType type : errorTypesNextGap()) {
		if (errorProbabilities[type] > bestProb) {
			bestProb = errorProbabilities[type];
			bestType = type;
		}
	}
	return std::make_pair(bestType, bestProb);
}

std::string kmerAfterError(const std::string &kmer, ErrorType error, size_t posOfError) {
	std::string res = kmer;
	if (error == ErrorType::SUB_FROM_A) {
		res[posOfError] = 'A';
	} else if (error == ErrorType::SUB_FROM_C) {
		res[posOfError] = 'C';
	} else if (error == ErrorType::SUB_FROM_G) {
		res[posOfError] = 'G';
	} else if (error == ErrorType::SUB_FROM_T) {
		res[posOfError] = 'T';
	} else if (error == ErrorType::INSERTION) {
		res = kmer.substr(0, posOfError) + kmer.substr(posOfError + 1, kmer.size());
	} else if (error == ErrorType::DEL_OF_A) {
		res = kmer.substr(0, posOfError + 1) + "A" + kmer.substr(posOfError + 1, kmer.size());
	} else if (error == ErrorType::DEL_OF_C) {
		res = kmer.substr(0, posOfError + 1) + "C" + kmer.substr(posOfError + 1, kmer.size());
	} else if (error == ErrorType::DEL_OF_G) {
		res = kmer.substr(0, posOfError + 1) + "G" + kmer.substr(posOfError + 1, kmer.size());
	} else if (error == ErrorType::DEL_OF_T) {
		res = kmer.substr(0, posOfError + 1) + "T" + kmer.substr(posOfError + 1, kmer.size());
	} else { // multideletion or chimeric break
		res = kmer.substr(0, posOfError + 1) + "_" + kmer.substr(posOfError + 1, kmer.size());
	}

	return res;
}

// TODO: FIXME: Improve handling of multideletions here
void correctKmer(const std::string &kmer, size_t kmerStartPos, CorrectedRead &corr, ErrorProfileUnit &errorProfile,
		KmerClassificationUnit &kmerClassifier) {
	std::vector<std::unordered_map<ErrorType, double> > probs = errorProfile.getReadErrorProbabilitiesPartial(
			corr.correctedRead, kmerStartPos, kmerStartPos + kmer.size());

	// sort the possible corrections based on their probability
	std::vector<std::vector<std::pair<ErrorType, double> > > probRankings;
	for (size_t i = 0; i < kmer.size(); ++i) {
		std::vector<std::pair<ErrorType, double> > rankings;
		for (auto kv : probs[i]) {
			rankings.push_back(std::pair<ErrorType, double>(kv.first, kv.second));
		}
		// sort rankings by descending kv.second
		std::sort(rankings.begin(), rankings.end(), [](const std::pair<ErrorType,double> &left, const std::pair<ErrorType,double> &right) {
		    return left.second > right.second;
		});
		probRankings.push_back(rankings);
	}

	for (size_t i = 0; i < kmer.size(); ++i) {
		for (size_t j = 0; j < probRankings[i].size(); ++j) {
			ErrorType bestError = probRankings[i][j].first;
			if (bestError == ErrorType::MULTIDEL) continue; // ignore Multidels for now... TODO: change this
			std::string correctedKmer = kmerAfterError(kmer, bestError, i);
			KmerType type = kmerClassifier.classifyKmer(correctedKmer);

			// extend the k-mer if it is repetitive now
			size_t inc = 1;
			while (type == KmerType::REPEAT
					&& kmerStartPos + correctedKmer.size() + inc + 1 < corr.correctedRead.sequence.size()) {
				correctedKmer += corr.correctedRead.sequence[kmerStartPos + correctedKmer.size() + inc];
				correctedKmer += corr.correctedRead.sequence[kmerStartPos + correctedKmer.size() + inc + 1];
				type = kmerClassifier.classifyKmer(correctedKmer);
				inc += 2;
			}

			// check if the k-mer is fine now
			if (type != KmerType::UNTRUSTED) {
				// found the correction. Apply the correction to the corrected read,
				corr.applyCorrection(bestError, i, probRankings[i][j].second);
				return;
			}
		}
	}
}

// TODO FIXME: This is currently not detecting multideletions. Also, the read is currently not correctly covered by k-mers.
CorrectedRead precorrectRead_KmerBased(const FASTQRead &fastqRead, ErrorProfileUnit &errorProfile,
		KmerClassificationUnit &kmerClassifier) {
	CorrectedRead corr(fastqRead);
	// cover the read with k-mers and correct them
	size_t i = 0;
	while (i < corr.correctedRead.sequence.size()) {
		std::string kmerString = corr.correctedRead.sequence.substr(i, kmerClassifier.getMinKmerSize());
		if (kmerString.find("_") != std::string::npos) { // kmer contains multidel, thus it should be ignored here
			continue;
		}
		KmerType kmerType = kmerClassifier.classifyKmer(kmerString);
		size_t j = 1;
		while (kmerType == KmerType::REPEAT && i + kmerString.size() + j + 1 < corr.correctedRead.sequence.size()) {
			if (corr.correctedRead.sequence[i + kmerString.size() + j] == '_'
					|| corr.correctedRead.sequence[i + kmerString.size() + j + 1] == '_') {
				break;
			}
			kmerString += corr.correctedRead.sequence[i + kmerString.size() + j];
			kmerString += corr.correctedRead.sequence[i + kmerString.size() + j + 1];
			kmerType = kmerClassifier.classifyKmer(kmerString);
		}
		if (kmerType == KmerType::UNTRUSTED) {
			correctKmer(kmerString, i, corr, errorProfile, kmerClassifier);
		}
		i++;
	}

	return corr;
}

CorrectedRead precorrectRead_Naive(const FASTQRead &fastqRead, ErrorProfileUnit &errorProfile) {
	CorrectedRead corr(fastqRead);
	int i = 0;
	while (i < (int) corr.correctedRead.sequence.size()) {
		auto probs = errorProfile.getErrorProbabilities(corr.correctedRead, i);
		std::pair<ErrorType, double> bestCurrent = mostLikelyCurrentBase(probs);
		if (bestCurrent.first != ErrorType::CORRECT) {

			if (!((corr.correctedRead.sequence[i] == 'A' && bestCurrent.first == ErrorType::SUB_FROM_A)
					|| (corr.correctedRead.sequence[i] == 'C' && bestCurrent.first == ErrorType::SUB_FROM_C)
					|| (corr.correctedRead.sequence[i] == 'G' && bestCurrent.first == ErrorType::SUB_FROM_G)
					|| (corr.correctedRead.sequence[i] == 'T' && bestCurrent.first == ErrorType::SUB_FROM_T))) {

				corr.applyCorrection(bestCurrent.first, i, bestCurrent.second);
				if (bestCurrent.first == ErrorType::INSERTION) {
					assert(i > 0); // there should be no insertion at the beginning of a read
					i--;
				}
			}
		}

		probs = errorProfile.getErrorProbabilities(corr.correctedRead, i);
		std::pair<ErrorType, double> bestNext = mostLikelyNextGap(probs);
		if (bestNext.first != ErrorType::NODEL) {
			corr.applyCorrection(bestNext.first, i, bestNext.second);
		}

		++i;
	}

	return corr;
}

CorrectedRead postcorrectRead_Multidel(const FASTQRead &fastqRead, ErrorProfileUnit &errorProfile,
		KmerClassificationUnit &kmerClassifier) {
	throw std::runtime_error("THIS IS NOT IMPLEMENTED YET.");
	CorrectedRead corr(fastqRead);
	int i = 1;
	while (i < (int) corr.correctedRead.sequence.size()) {
		if (corr.correctedRead.sequence[i] == '_') { // encountered a multideletion... trying to fix it.
			auto probs = errorProfile.getErrorProbabilities(corr.correctedRead, i - 1);
			// TODO: continue implementation
		}

		/*std::pair<ErrorType, double> bestCurrent = mostLikelyCurrentBase(probs);
		 if (bestCurrent.first != ErrorType::CORRECT) {
		 corr.applyCorrection(bestCurrent.first, i, bestCurrent.second);
		 if (bestCurrent.first == ErrorType::INSERTION) {
		 assert(i > 0); // there should be no insertion at the beginning of a read
		 i--;
		 }
		 }

		 probs = errorProfile->getErrorProbabilities(corr.correctedRead, i);
		 std::pair<ErrorType, double> bestNext = mostLikelyNextGap(probs);
		 if (bestNext.first != ErrorType::NODEL) {
		 corr.applyCorrection(bestNext.first, i, bestNext.second);
		 }*/

		++i;
	}

	return corr;
}
