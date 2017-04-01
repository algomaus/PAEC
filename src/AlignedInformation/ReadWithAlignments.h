/*
 * ReadWithAlignments.h
 *
 *  Created on: Jan 31, 2017
 *      Author: sarah
 */

#pragma once

#include <seqan/bam_io.h>

#include <vector>

#include "CorrectedReadAligned.h"

/*
 * TODO: Bei einem Chimeric break dürfte die beginPos des Reads unterschiedlich sein!!!
 */

class ReadWithAlignments {
public:
	ReadWithAlignments(std::vector<seqan::BamAlignmentRecord> &bamRecords);
	CorrectedReadAligned retrieveCorrectedRead(const seqan::Dna5String &referenceGenome);

	seqan::Dna5String myReverseComplement(const seqan::Dna5String &read);

	seqan::Dna5 reverseComplementBase(const seqan::Dna5 &base);

	seqan::CharString name;seqan::Dna5String seq;seqan::CharString qual;
	size_t beginPos;

	std::vector<seqan::BamAlignmentRecord> records;
private:
	void extractErrors(const seqan::Dna5String &referenceGenome);
	CorrectedReadAligned correctedRead;
};
