/*
 * GenomeReader.cpp
 *
 *  Created on: Jan 31, 2017
 *      Author: sarah
 */

#include "GenomeReader.h"

#include <iostream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

GenomeReader::GenomeReader() {}

std::string GenomeReader::readGenome(const std::string &genomeFilename) {
	seqan::Dna5String referenceSequence;
	seqan::CharString id;
	seqan::SeqFileIn seqFileIn(genomeFilename.c_str());
	seqan::readRecord(id, referenceSequence, seqFileIn);

	seqan::CharString ref = referenceSequence;
	return std::string(toCString(ref));
}
