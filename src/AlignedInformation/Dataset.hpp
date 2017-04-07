#pragma once

#include <fstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <limits>

#include <seqan/sequence.h>

#include "GenomeReader.h"
#include "../CoverageBias/GenomeType.h"

static const int MIN_KMER_SIZE = 15; // TODO: Does it make sense to introduce a minimum k-mer size??? Yes it does. For example, I use it in FeatureExtractorCurrentBase.hpp.
static const int MAX_KMER_SIZE = 35;

class Dataset {
public:
	Dataset() {
		numReads = 0;
		genomeType = GenomeType::CIRCULAR;
		genomeSize = 0;
		minReadLength = 0;
		hasQualityScores = false;
		maxReadLength = 0;
	}

	Dataset(const std::string &sra) {
		hasQualityScores = true;
		if (sra == "SRR396537") {
			readsOnlyFileName =
					"data/e_coli_k12_mg1655/Illumina_Genome_Analyzer_IIx/Nextera/GA091120/SRR396537/SRR396537_without_adapters.fastq.readsOnly.txt";
			readsFileName =
					"data/e_coli_k12_mg1655/Illumina_Genome_Analyzer_IIx/Nextera/GA091120/SRR396537/SRR396537_without_adapters.fastq";
			referenceFileName = "data/e_coli_k12_mg1655/reference.fasta";
			readAlignmentsFileName =
					"data/e_coli_k12_mg1655/Illumina_Genome_Analyzer_IIx/Nextera/GA091120/SRR396537/SRR396537_without_adapters.sorted.bam";
			plotPath = "plots/Illumina/SRR396537/SRR396537_";
			genomeType = GenomeType::CIRCULAR;
		} else if (sra == "SRR396536") {
			readsOnlyFileName =
					"data/e_coli_k12_mg1655/Illumina_Genome_Analyzer_IIx/Nextera/GA091120/SRR396536/SRR396536_without_adapters.fastq.readsOnly.txt";
			readsFileName =
					"data/e_coli_k12_mg1655/Illumina_Genome_Analyzer_IIx/Nextera/GA091120/SRR396536/SRR396536_without_adapters.fastq";
			referenceFileName = "data/e_coli_k12_mg1655/reference.fasta";
			readAlignmentsFileName =
					"data/e_coli_k12_mg1655/Illumina_Genome_Analyzer_IIx/Nextera/GA091120/SRR396536/SRR396536_without_adapters.sorted.bam";
			plotPath = "plots/Illumina/SRR396536/SRR396536_";
			genomeType = GenomeType::CIRCULAR;
		} else if (sra == "SRR1284073") {
			readsOnlyFileName = "data/e_coli_k12_mg1655/PacBio/SRR1284073.fastq.readsOnly.txt";
			readsFileName =
					"data/e_coli_k12_mg1655/PacBio/SRR1284073.fastq";
			referenceFileName = "data/e_coli_k12_mg1655/reference.fasta";
			readAlignmentsFileName =
					"data/e_coli_k12_mg1655/PacBio/SRR1284073.sorted.bam";
			plotPath = "plots/PacBio/SRR1284073/SRR1284073_";
			genomeType = GenomeType::CIRCULAR;
		} else if (sra == "ebola_pacbio_simulated") {
			readsOnlyFileName = "data/Simulated Datasets/Ebola/PacBio/ebola_pacbio_simulated.fastq.readsOnly.txt";
			readsFileName = "data/Simulated Datasets/Ebola/PacBio/ebola_pacbio_simulated.fastq";
			referenceFileName = "data/Simulated Datasets/Ebola/reference.fasta";
			readAlignmentsFileName = "data/Simulated Datasets/Ebola/PacBio/ebola_pacbio_simulated.bam";
			plotPath = "plots/PacBio/ebola_simulated/ebola_pacbio_simulated_";
			genomeType = GenomeType::LINEAR;
		} else if (sra == "ebola_illumina_simulated") {
			readsOnlyFileName = "data/Simulated Datasets/Ebola/Illumina/ebola_illumina_simulated.fastq.readsOnly.txt";
			readsFileName = "data/Simulated Datasets/Ebola/Illumina/ebola_illumina_simulated.fastq";
			referenceFileName = "data/Simulated Datasets/Ebola/reference.fasta";
			readAlignmentsFileName = "data/Simulated Datasets/Ebola/Illumina/ebola_illumina_simulated.bam";
			plotPath = "plots/Illumina/ebola_simulated/ebola_illumina_simulated_";
			genomeType = GenomeType::LINEAR;
		} else if (sra == "ecoli_pacbio_simulated") {
			readsOnlyFileName = "data/Simulated Datasets/Ecoli/PacBio/ecoli_pacbio_simulated.fastq.readsOnly.txt";
			readsFileName = "data/Simulated Datasets/Ecoli/PacBio/ecoli_pacbio_simulated.fastq";
			referenceFileName = "data/Simulated Datasets/Ecoli/reference.fasta";
			readAlignmentsFileName = "data/Simulated Datasets/Ecoli/PacBio/ecoli_pacbio_simulated.bam";
			plotPath = "plots/PacBio/ecoli_simulated/ecoli_pacbio_simulated_";
			genomeType = GenomeType::CIRCULAR;
		} else if (sra == "ecoli_illumina_simulated") {
			readsOnlyFileName = "data/Simulated Datasets/Ecoli/Illumina/ecoli_illumina_simulated.fastq.readsOnly.txt";
			readsFileName = "data/Simulated Datasets/Ecoli/Illumina/ecoli_illumina_simulated.fastq";
			referenceFileName = "data/Simulated Datasets/Ecoli/reference.fasta";
			readAlignmentsFileName = "data/Simulated Datasets/Ecoli/Illumina/ecoli_illumina_simulated.bam";
			plotPath = "plots/Illumina/ecoli_simulated/ecoli_illumina_simulated_";
			genomeType = GenomeType::CIRCULAR;
		} else {
			throw "Unknown dataset!";
		}
		name = sra;
		GenomeReader reader;
		genome = reader.readGenome(referenceFileName);
		genomeSize = length(genome);
		countReadLengths();
	}
	std::string readsOnlyFileName;
	std::string referenceFileName;
	std::string readsFileName;
	std::string readAlignmentsFileName;
	std::string plotPath;
	seqan::Dna5String genome;
	size_t genomeSize;
	GenomeType genomeType;
	std::unordered_map<size_t, size_t> readLengths;
	std::shared_ptr<std::unordered_map<size_t, size_t> > readLengthsPtr;
	size_t minReadLength;
	size_t maxReadLength;
	size_t numReads;
	bool hasQualityScores;
	double acceptProb = 1.0;
	std::string name;
private:
	void countReadLengths() {
		numReads = 0;
		minReadLength = std::numeric_limits<size_t>::infinity();
		maxReadLength = 0;
		std::ifstream infile(readsOnlyFileName);
		std::string line;
		while (std::getline(infile, line)) {
			size_t s = line.size();
			if (s > 0) {
				if (readLengths.find(s) == readLengths.end()) {
					readLengths[s] = 1;
				} else {
					readLengths[s]++;
				}
				minReadLength = std::min(minReadLength, s);
				maxReadLength = std::max(maxReadLength, s);
				numReads++;
			}
		}
		readLengthsPtr = std::make_shared<std::unordered_map<size_t, size_t> >(readLengths);
	}
};
