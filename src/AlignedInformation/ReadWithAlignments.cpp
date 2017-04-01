/*
 * ReadWithAlignments.cpp
 *
 *  Created on: Jan 31, 2017
 *      Author: sarah
 */

#include "ReadWithAlignments.h"

#include <cassert>
#include <iostream>
#include <string>
#include <algorithm>

#include "../ErrorType.h"
#include "mdTag.hpp"

CorrectedReadAligned ReadWithAlignments::retrieveCorrectedRead(const seqan::Dna5String &referenceGenome) {
	extractErrors(referenceGenome);
	return correctedRead;
}

ReadWithAlignments::ReadWithAlignments(std::vector< seqan::BamAlignmentRecord> &bamRecords) {
	records = bamRecords;

	// go through records, try to find one which has no hard-clipping
	bool found = false;
	for (seqan::BamAlignmentRecord record : records) {
		found = true;
		seqan::String<seqan::CigarElement<char> > cigar = record.cigar;
		for (size_t i = 0; i < length(cigar); ++i) {
			if (cigar[i].operation == 'H') {
				found = false;
				break;
			}
		}
		if (found) { // found a record which has no hard clipping
			FASTQRead fastQread;

			fastQread.id = toCString(record.qName);
			seq = "";
			qual = "";

			// copy the read
			for (size_t i = 0; i < length(record.seq); ++i) {
				seq += record.seq[i];
			}

			// if the read is reverse-complemented, reverse-complement it again to be normal.
			if (hasFlagRC(record)) {
				seq = myReverseComplement(seq);
				for (int i = length(record.qual) - 1; i >= 0; --i) {
					qual += record.qual[i];
				}
			} else {
				for (size_t i = 0; i < length(record.qual); ++i) {
					qual += record.qual[i];
				}
			}

			seqan::CharString sequenceAsCharString = seq;
			fastQread.sequence = toCString(sequenceAsCharString);
			fastQread.quality = toCString(qual);
			correctedRead = CorrectedReadAligned(fastQread, record.beginPos);
			break;
		}
	}
	if (!found) {
		throw "ERROR! Read only occurs with hard clipping.\n";
	}
}

seqan::Dna5String ReadWithAlignments::myReverseComplement(const seqan::Dna5String &read) {
	seqan::Dna5String rcRead = "";
	for (int i = length(read) - 1; i >= 0; --i) {
		if (read[i] == seqan::Dna5('A')) {
			rcRead += seqan::Dna5('T');
		} else if (read[i] == seqan::Dna5('C')) {
			rcRead += seqan::Dna5('G');
		} else if (read[i] == seqan::Dna5('G')) {
			rcRead += seqan::Dna5('C');
		} else if (read[i] == seqan::Dna5('T')) {
			rcRead += seqan::Dna5('A');
		} else if (read[i] == seqan::Dna5('N')) {
			rcRead += seqan::Dna5('N');
		} else {
			std::cout << "Strange base: " << read[i] << std::endl;
			throw "Error. Base cannot be reverse complemented.";
		}
	}
	return rcRead;
}

void ReadWithAlignments::extractErrors(const seqan::Dna5String &genome) {
	// ignore searching errors in unmapped reads
	if (hasFlagUnmapped(records[0])) {
		return;
	}

	// ignore searching errors in soft-clipped non-chimeric reads
	if (records.size() == 1) {
		for (size_t i = 0; i < length(records[0].cigar); ++i) {
			if (records[0].cigar[i].operation == 'S') {
				return;
			}
		}
	}

	// ignore searching errors in reads with chimeric breaks
	/*if (records.size() > 1) {
	 return;
	 }*/

	std::vector<CorrectionAligned> corrections;

	for (seqan::BamAlignmentRecord record : records) {
		// Extract MD tag, CIGAR string, total read length (including unmapped regions)
		MDTag mdTag;
		extractMDTag(record, mdTag);
		seqan::String< seqan::CigarElement<char> > cigar = record.cigar;
		unsigned readLength = 0;

		for (unsigned i = 0; i < length(cigar); ++i) {
			if (cigar[i].operation != 'D') {
				readLength += cigar[i].count;
			}
		}
		assert(readLength == length(alignedRead.read));

		unsigned positionInRead = 0;
		unsigned hardClippedBases = 0;
		unsigned softClippedBases = 0;
		unsigned insertedBases = 0;
		unsigned deletedBases = 0;
		unsigned actMDindex = 0;

		for (unsigned i = 0; i < length(cigar); ++i) {
			if (cigar[i].operation == 'H') { // hard clipping
				unsigned nucleotidePositionInReference = record.beginPos + positionInRead - insertedBases
						- softClippedBases + deletedBases;
				unsigned realPositionInRead = positionInRead + hardClippedBases; // the position of the first base of the chimeric break

				if (realPositionInRead == 0) { // in this case, take the position of the last base of the chimeric break
					realPositionInRead += cigar[i].count;
				}

				// TODO: What if multiple chimeric breaks happen in a read? ... It should be fine, too. As long as there aren't multiple hard clipped regions in a single record.

				hardClippedBases += cigar[i].count;
				if (hasFlagRC(record)) {
					realPositionInRead = readLength - realPositionInRead - 1;
				}

				corrections.push_back(
						CorrectionAligned(nucleotidePositionInReference,
								Correction(realPositionInRead, "", "", ErrorType::MULTIDEL))); // because we treat chimeric breaks as deletion of multiple bases
			} else if (cigar[i].operation == 'S') { // soft clipping
				positionInRead += cigar[i].count;
				softClippedBases += cigar[i].count;
			} else if (cigar[i].operation == 'M') { // match ... check MD tag for substitutions
				for (size_t j = actMDindex; j < mdTag.size(); ++j) {
					// (because insertions are ignored in the MD tag)
					size_t cleanedMDPosition = mdTag[j].position + insertedBases + softClippedBases;
					if (cleanedMDPosition >= positionInRead && cleanedMDPosition <= positionInRead + cigar[i].count) {
						// Find nucleotide in reference genome that has been substituted
						unsigned nucleotidePositionInRead = cleanedMDPosition;
						seqan::Dna5 nucleotideInRead = record.seq[nucleotidePositionInRead];
						unsigned nucleotidePositionInReference = record.beginPos + nucleotidePositionInRead
								- insertedBases - softClippedBases + deletedBases;
						seqan::Dna5 nucleotideInReference = genome[nucleotidePositionInReference];
						unsigned realPositionInRead = nucleotidePositionInRead + hardClippedBases;
						if (hasFlagRC(record)) {
							realPositionInRead = readLength - realPositionInRead - 1;
							nucleotideInRead = reverseComplementBase(nucleotideInRead);
							nucleotideInReference = reverseComplementBase(nucleotideInReference);
						}
						// TODO: Ist das die korrekte Position im Read?
						assert(nucleotideInRead == alignedRead.read[realPositionInRead]);

						std::string fromBase = "";
						fromBase += nucleotideInRead;
						if (nucleotideInReference == 'A') {
							assert(fromBase[0] != 'A');
							corrections.push_back(
									CorrectionAligned(nucleotidePositionInReference,
											Correction(realPositionInRead, fromBase, "A", ErrorType::SUB_FROM_A)));
						} else if (nucleotideInReference == 'C') {
							assert(fromBase[0] != 'C');
							corrections.push_back(
									CorrectionAligned(nucleotidePositionInReference,
											Correction(realPositionInRead, fromBase, "C", ErrorType::SUB_FROM_C)));
						} else if (nucleotideInReference == 'G') {
							assert(fromBase[0] != 'G');
							corrections.push_back(
									CorrectionAligned(nucleotidePositionInReference,
											Correction(realPositionInRead, fromBase, "G", ErrorType::SUB_FROM_G)));
						} else if (nucleotideInReference == 'T') {
							assert(fromBase[0] != 'T');
							corrections.push_back(
									CorrectionAligned(nucleotidePositionInReference,
											Correction(realPositionInRead, fromBase, "T", ErrorType::SUB_FROM_T)));
						} else {
							throw "N in reference genome!";
						}

						actMDindex++;
					} else {
						break;
					}
				}
				positionInRead += cigar[i].count;
			} else if (cigar[i].operation == 'I') { // insertion
				for (size_t j = 0; j < cigar[i].count; ++j) {
					unsigned nucleotidePositionRead = positionInRead;
					seqan::Dna5 nucleotideInRead = record.seq[nucleotidePositionRead];
					unsigned nucleotidePositionInReference = record.beginPos + nucleotidePositionRead - insertedBases
							- softClippedBases + deletedBases;
					size_t realPositionInRead = positionInRead + hardClippedBases;

					if (hasFlagRC(record)) {
						realPositionInRead = readLength - realPositionInRead - 1;
						nucleotideInRead = reverseComplementBase(nucleotideInRead);
					}
					std::string fromBase = "";
					fromBase += nucleotideInRead;
					corrections.push_back(
							CorrectionAligned(nucleotidePositionInReference,
									Correction(realPositionInRead, fromBase, "", ErrorType::INSERTION)));
					positionInRead++; // TODO: Davor oder danach???
				}
				insertedBases += cigar[i].count;

			} else if (cigar[i].operation == 'D') { // deletion
				unsigned nucleotidePositionRead = positionInRead;

				size_t realPositionInRead = positionInRead + hardClippedBases;
				if (hasFlagRC(record)) {
					realPositionInRead = readLength - realPositionInRead - 1;
				}
				std::string fromBase = "";
				fromBase += correctedRead.correctedRead.sequence[realPositionInRead];
				std::string toBases = "";
				unsigned nucleotidePositionInReference = record.beginPos + nucleotidePositionRead - insertedBases
						- softClippedBases + deletedBases;
				seqan::Dna5 nucleotideInReference = genome[nucleotidePositionInReference];

				for (size_t j = 0; j < cigar[i].count; ++j) {
					deletedBases++;
					nucleotideInReference = genome[nucleotidePositionInReference + j + 1];
					if (hasFlagRC(record)) {
						nucleotideInReference = reverseComplementBase(nucleotideInReference);
					}
					toBases += nucleotideInReference;
				}
				assert(toBases.size() == cigar[i].count);

				if (cigar[i].count == 1) {
					if (nucleotideInReference == 'A') {
						corrections.push_back(
								CorrectionAligned(nucleotidePositionInReference,
										Correction(realPositionInRead, fromBase, fromBase + "A", ErrorType::DEL_OF_A)));
					} else if (nucleotideInReference == 'C') {
						corrections.push_back(
								CorrectionAligned(nucleotidePositionInReference,
										Correction(realPositionInRead, fromBase, fromBase + "C", ErrorType::DEL_OF_C)));
					} else if (nucleotideInReference == 'G') {
						corrections.push_back(
								CorrectionAligned(nucleotidePositionInReference,
										Correction(realPositionInRead, fromBase, fromBase + "G", ErrorType::DEL_OF_G)));
					} else if (nucleotideInReference == 'T') {
						corrections.push_back(
								CorrectionAligned(nucleotidePositionInReference,
										Correction(realPositionInRead, fromBase, fromBase + "T", ErrorType::DEL_OF_T)));
					} else {
						throw "N in reference genome!";
					}
				} else {
					corrections.push_back(
							CorrectionAligned(nucleotidePositionInReference,
									Correction(realPositionInRead, fromBase, fromBase + toBases, ErrorType::MULTIDEL)));
				}

			} else if (cigar[i].operation == 'P' || cigar[i].operation == 'N') {
				std::cout << "ERROR! There are letters I don't understand yet!" << cigar[i].operation << std::endl;
			}
		}
	}

	std::sort(corrections.begin(), corrections.end(), [](const CorrectionAligned& a, const CorrectionAligned& b)
	{
		return a.correction.positionInRead < b.correction.positionInRead;
	});
	for (size_t i = 0; i < corrections.size(); ++i) {
		correctedRead.applyCorrection(corrections[i]);
	}
}

seqan::Dna5 ReadWithAlignments::reverseComplementBase(const seqan::Dna5 &base) {
	if (base == seqan::Dna5('A')) {
		return seqan::Dna5('T');
	} else if (base == seqan::Dna5('C')) {
		return seqan::Dna5('G');
	} else if (base == seqan::Dna5('G')) {
		return seqan::Dna5('C');
	} else if (base == seqan::Dna5('T')) {
		return seqan::Dna5('A');
	} else {
		return seqan::Dna5('N');
	}
}
