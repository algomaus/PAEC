#include "CoverageBiasUnit.h"

#include "../AlignedInformation/BAMIterator.h"
#include "../AlignedInformation/ReadWithAlignments.h"
#include "../FASTQIterator.h"
#include "../FASTQRead.h"

#include "../external/IntervalTree.h"
#include "../Plotter.hpp"

// TODO: Adapt the count in the reference genome to also deal with the special case of a circular genome

#define MIN_TRUSTED_COUNT 2

CoverageBiasUnit::CoverageBiasUnit() {
	minKmerSize = 1;
	genomeSize = 1;
	gc_step = 0.1;
	covBiasType = CoverageBiasType::IGNORE_BIAS;
}

CoverageBiasUnit::CoverageBiasUnit(CoverageBiasType coverageBiasType, size_t estimatedGenomeSize,
		PerfectUniformSequencingModel &pusm_ref) {
	minKmerSize = std::floor(log(200 * estimatedGenomeSize) / log(4));
	if (minKmerSize % 2 == 0) {
		minKmerSize++;
	}
	std::cout << "minKmerSize = " << minKmerSize << "\n";
	genomeSize = estimatedGenomeSize;

	gc_step = 1.0 / minKmerSize;
	assert(minKmerSize % 2 == 1); // ensure that reverse complement of a kmer is never the same as the kmer
	biases.resize(1.0 / gc_step + 1);

	pusm = std::make_shared<PerfectUniformSequencingModel>(pusm_ref);
	covBiasType = coverageBiasType;
}

void CoverageBiasUnit::loadBias(const std::string &filepath) {
	std::ifstream infile(filepath, std::ios::binary);
	if (!infile.good()) {
		throw std::runtime_error("The file " + filepath + " does not exist!");
	}
	cereal::BinaryInputArchive iarchive(infile);
	CoverageBiasUnit mcb;
	iarchive(mcb);
	minKmerSize = mcb.minKmerSize;
	gc_step = mcb.gc_step;
	biases = mcb.biases;
}

void CoverageBiasUnit::storeBias(const std::string &filepath) {
	std::ofstream outfile(filepath, std::ios::binary);
	cereal::BinaryOutputArchive oarchive(outfile);
	oarchive(*this);
}

double CoverageBiasUnit::getBias(const std::string &kmer) {
	if (kmer.empty()) {
		throw std::runtime_error("The kmer is empty!");
	}
	if (covBiasType == CoverageBiasType::IGNORE_BIAS) {
		return 1.0;
	} else {
		double gc = computeGCContent(kmer);
		size_t idxMin = std::floor(gc / gc_step);
		size_t idxMax = idxMin + 1;
		double gcMin = idxMin * gc_step;
		double gcMax = idxMax * gc_step;
		double biasMin = biases[idxMin];
		double biasMax = biases[idxMax];
		// linear interpolation, see https://en.wikipedia.org/wiki/Interpolation
		return biasMin + (biasMax - biasMin) / (gcMax - gcMin) * (gc - gcMin);
	}
}

void CoverageBiasUnit::learnBiasFromReferenceAlignment(const seqan::Dna5String &referenceGenome,
		KmerCounter &referenceCounter, const std::string &alignmentsFilename,
		const std::unordered_map<size_t, size_t> &readLengths) {
	if (covBiasType == CoverageBiasType::IGNORE_BIAS)
		return;
	if (covBiasType != CoverageBiasType::MEDIAN_BIAS_ALIGNMENT) {
		throw std::runtime_error("Wrong coverage bias type");
	}

	std::cout << "Constructing interval counts...\n";

	//auto intervalCounts = std::vector<double>( static_cast<size_t>( length(referenceGenome) ), 0.0);

//	std::vector<double> intervalCounts;
//	for (size_t i = 0; i < length(referenceGenome); ++i) {
//		intervalCounts.push_back(0);
//	}

	std::vector<Interval<double> > intervals;

	BAMIterator it(alignmentsFilename);

	while (it.hasReadsLeft()) {
		std::vector<ReadWithAlignments> readsBuffer = it.next(100);
		for (ReadWithAlignments read : readsBuffer) {
			for (seqan::BamAlignmentRecord record : read.records) {
				if (!hasFlagUnmapped(record)) {
					size_t intervalStart = record.beginPos;
					size_t intervalEnd = intervalStart + length(record.seq) - 1;
					//std::cout << "intervalStart: " << intervalStart << ", intervalEnd: " << intervalEnd << "\n";
					assert(intervalStart <= intervalEnd);
					assert(intervalEnd <= length(referenceGenome));
					intervals.push_back(Interval<double>(intervalStart, intervalEnd, 1.0 / read.records.size()));
				}
			}
		}
	}

	std::cout << "intervals.size() = " << intervals.size() << "\n";

	IntervalTree<double> tree;
	tree = IntervalTree<double>(intervals);

	std::vector<std::vector<double> > allBiases;
	allBiases.resize(biases.size());

	std::cout << "Learning coverage biases from reference genome and read dataset, using interval tree...\n";
	size_t min_progress = 0;
	std::string kmerString = "_";
	size_t numGC = 0;
	for (size_t i = 0; i < minKmerSize - 1; ++i) {
		kmerString += referenceGenome[i];
		if (referenceGenome[i] == 'G' || referenceGenome[i] == 'C') {
			numGC++;
		}
	}

	for (size_t i = minKmerSize - 1; i < length(referenceGenome); ++i) {
		if (kmerString[0] == 'G' || kmerString[0] == 'C') {
			numGC--;
		}
		kmerString = kmerString.substr(1, kmerString.size() - 1); // delete first character
		kmerString += referenceGenome[i]; // append current character
		if (referenceGenome[i] == 'G' || referenceGenome[i] == 'C') {
			numGC++;
		}

		size_t occRef = referenceCounter.countKmerNoRC(kmerString);

		double gc = (double) numGC / minKmerSize;
		size_t gcIndex = gc / gc_step;

		double countObserved = 0.0;

		auto coveringIntervals = tree.findCovering((i + 1) - minKmerSize, i);
		for (size_t i = 0; i < coveringIntervals.size(); ++i) {
			countObserved += coveringIntervals[i].value;
		}

		if (countObserved > MIN_TRUSTED_COUNT) { // if this condition is left out, the PacBio dataset will get a median coverage bias of 0.
			double countExpected = pusm->expectedCount(kmerString).first;
			countExpected *= occRef;
			double bias = countObserved / countExpected;
			allBiases[gcIndex].push_back(bias);
		}

		double progress = 100.0 * (double) i / (length(referenceGenome) - minKmerSize + 1);
		if (progress >= min_progress) {
			std::cout << progress << "%\n";
			min_progress++;
		}
	}

	for (size_t i = 0; i < allBiases.size(); ++i) {
		std::sort(allBiases[i].begin(), allBiases[i].end());
		size_t sizeHalved = allBiases[i].size() / 2;
		if (allBiases[i].size() == 0) {
			biases[i] = 0.0;
		} else {
			if (allBiases[i].size() % 2 == 1) {
				biases[i] = allBiases[i][sizeHalved + 1];
			} else {
				biases[i] = (allBiases[i][sizeHalved] + allBiases[i][sizeHalved + 1]) / 2;
			}
		}
	}
	fixEmptyBiases();

	std::cout << "Finished learning coverage biases from reference genome and read dataset, using interval tree.\n";
}

void CoverageBiasUnit::learnBiasFromReferenceMatches(const seqan::Dna5String &referenceGenome,
		KmerCounter &referenceCounter, KmerCounter &readsCounter,
		const std::unordered_map<size_t, size_t> &readLengths) {
	if (covBiasType == CoverageBiasType::IGNORE_BIAS)
		return;
	if (covBiasType != CoverageBiasType::MEDIAN_BIAS_REFERENCE) {
		throw std::runtime_error("Wrong coverage bias type");
	}
	std::vector<std::vector<double> > allBiases;
	allBiases.resize(biases.size());

	std::cout << "Learning coverage biases from reference genome and read dataset, exact matches only...\n";
	size_t min_progress = 0;
	std::string kmerString = "_";
	size_t numGC = 0;
	for (size_t i = 0; i < minKmerSize - 1; ++i) {
		kmerString += referenceGenome[i];
		if (referenceGenome[i] == 'G' || referenceGenome[i] == 'C') {
			numGC++;
		}
	}

	for (size_t i = minKmerSize - 1; i < length(referenceGenome); ++i) {
		if (kmerString[0] == 'G' || kmerString[0] == 'C') {
			numGC--;
		}
		kmerString = kmerString.substr(1, kmerString.size() - 1); // delete first character
		kmerString += referenceGenome[i]; // append current character
		if (referenceGenome[i] == 'G' || referenceGenome[i] == 'C') {
			numGC++;
		}

		size_t occRef = referenceCounter.countKmerNoRC(kmerString);
		double gc = (double) numGC / minKmerSize;
		size_t gcIndex = gc / gc_step;
		size_t countObserved = readsCounter.countKmer(kmerString);
		if (countObserved > MIN_TRUSTED_COUNT) { // if this condition is left out, the PacBio dataset will get a median coverage bias of 0.
			double countExpected = pusm->expectedCount(kmerString).first;
			countExpected *= occRef;
			double bias = (double) countObserved / countExpected;
			allBiases[gcIndex].push_back(bias);
		}

		double progress = 100.0 * (double) i / (length(referenceGenome) - minKmerSize + 1);
		if (progress >= min_progress) {
			std::cout << progress << "%\n";
			min_progress++;
		}
	}

	//#pragma omp parallel for
	for (size_t i = 0; i < allBiases.size(); ++i) {
		std::sort(allBiases[i].begin(), allBiases[i].end());
		size_t sizeHalved = allBiases[i].size() / 2;
		if (allBiases[i].size() == 0) {
			biases[i] = 0; // TODO: Should this be another value?
		} else {
			if (allBiases[i].size() % 2 == 1) {
				biases[i] = allBiases[i][sizeHalved + 1];
			} else {
				biases[i] = (allBiases[i][sizeHalved] + allBiases[i][sizeHalved + 1]) / 2;
			}
		}
	}

	fixEmptyBiases();

	std::cout << "Finished learning coverage biases from reference genome and read dataset.\n";
}

void CoverageBiasUnit::learnBiasFromReadsOnly(const std::string &readsFilePath,
		KmerCounter &readsCounter, const std::unordered_map<size_t, size_t> &readLengths) {
	// go through all reads, go through all k-mers in those reads, for all unvisited k-mers update the coverage bias stuff
	if (covBiasType == CoverageBiasType::IGNORE_BIAS)
		return;
	if (covBiasType != CoverageBiasType::MEDIAN_BIAS_READS_ONLY) {
		throw std::runtime_error("Wrong coverage bias type");
	}
	std::vector<std::vector<double> > allBiases;
	allBiases.resize(biases.size());

	std::cout << "Learning coverage biases from read dataset only, exact matches only...\n";
	size_t min_progress = 0;

	std::unordered_set<std::string> visitedKmers;
	FASTQIterator it(readsFilePath);
	while (it.hasReadsLeft()) {
		std::vector<FASTQRead> reads = it.next(100);
		for (FASTQRead read : reads) {
			std::string seq = read.sequence;
			for (size_t i = 0; i < seq.size() - minKmerSize; ++i) {
				std::string kmer = read.sequence.substr(i, minKmerSize);
				if (visitedKmers.find(kmer) == visitedKmers.end()) {

					double gc = 0;
					for (size_t j = 0; j < kmer.size(); ++i) {
						if (kmer[i] == 'G' || kmer[i] == 'C') {
							gc++;
						}
					}
					gc = gc / kmer.size();
					size_t gcIndex = gc / gc_step;
					size_t countObserved = readsCounter.countKmer(kmer);
					if (countObserved > MIN_TRUSTED_COUNT) { // if this condition is left out, the PacBio dataset will get a median coverage bias of 0.
						double countExpected = pusm->expectedCount(kmer).first;
						double bias = (double) countObserved / countExpected;
						allBiases[gcIndex].push_back(bias);
					}

					visitedKmers.insert(kmer);
				}
			}
		}

		double progress = it.progress();
		if (progress >= min_progress) {
			std::cout << progress << "%\n";
			min_progress++;
		}
	}

	for (size_t i = 0; i < allBiases.size(); ++i) {
		std::sort(allBiases[i].begin(), allBiases[i].end());
		size_t sizeHalved = allBiases[i].size() / 2;
		if (allBiases[i].size() == 0) {
			biases[i] = 0; // TODO: Should this be another value?
		} else {
			if (allBiases[i].size() % 2 == 1) {
				biases[i] = allBiases[i][sizeHalved + 1];
			} else {
				biases[i] = (allBiases[i][sizeHalved] + allBiases[i][sizeHalved + 1]) / 2;
			}
		}
	}

	fixEmptyBiases();
}

void CoverageBiasUnit::fixEmptyBiases() {
	bool workLeft = true;
	while (workLeft) {
		workLeft = false;
		if (biases[biases.size() - 1] == 0.0) {
			if (biases.size() > 1) {
				bool foundNonZero = false;
				for (int i = biases.size() - 2; i >= 0; --i) {
					if (biases[i] != 0.0) {
						foundNonZero = true;
						biases[biases.size() - 1] = biases[i];
					}
				}
				if (!foundNonZero) {
					for (size_t i = 0; i < biases.size(); ++i) {
						biases[i] = 1.0;
					}
					std::cout << "All biases were zero. Fixed them by assuming no bias.\n";
					break;
				}
			} else {
				biases[0] = 1.0;
			}
			std::cout << "biases[" << (biases.size() - 1) * gc_step << "] was empty and got fixed.\n";
			workLeft = true;
		}

		if (biases[0] == 0.0) {
			if (biases.size() > 1) {
				bool foundNonZero = false;
				for (size_t i = 1; i < biases.size(); ++i) {
					if (biases[i] != 0.0) {
						foundNonZero = true;
						biases[0] = biases[i];
					}
				}
				if (!foundNonZero) {
					for (size_t i = 0; i < biases.size(); ++i) {
						biases[i] = 1.0;
					}
					std::cout << "All biases were zero. Fixed them by assuming no bias.\n";
					break;
				}
			} else {
				biases[0] = 1.0;
			}
			std::cout << "biases[0] was empty and got fixed.\n";
			workLeft = true;
		}

		for (size_t i = 1; i < biases.size() - 1; ++i) {
			if (biases[i] == 0.0) {
				biases[i] = (biases[i - 1] + biases[i + 1]) / 2;
				std::cout << "biases[" << i * gc_step << "] was empty and got fixed.\n";
				workLeft = true;
			}
		}
	}
}

void CoverageBiasUnit::printBias() { // TODO: Plot them using GNUplot...
	for (size_t i = 0; i < biases.size(); ++i) {
		std::cout << "biases[" << i * gc_step << "] = " << biases[i] << "\n";
	}
}

void CoverageBiasUnit::plotBias(const std::string &filename) {
	std::unordered_map<double, double> data;
	for (size_t i = 0; i < biases.size(); ++i) {
		data[gc_step * i] = biases[i];
	}
	plot(data, "G/C content", "coverage bias", filename);
}

double CoverageBiasUnit::computeGCContent(const std::string &sequence) {
	size_t gc = 0;
	for (size_t i = 0; i < sequence.size(); ++i) {
		if (sequence[i] == 'G' || sequence[i] == 'C') {
			gc++;
		}
	}
	return (double) gc / sequence.size();
}
