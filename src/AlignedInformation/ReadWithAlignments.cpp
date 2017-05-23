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
			fastQread.id = "@" + fastQread.id;
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

std::string reverseComplementString(const std::string &read) {
	std::string rcRead = "";
	for (int i = read.size() - 1; i >= 0; --i) {
		if (read[i] == 'A') {
			rcRead += 'T';
		} else if (read[i] == 'C') {
			rcRead += 'G';
		} else if (read[i] == 'G') {
			rcRead += 'C';
		} else if (read[i] == 'T') {
			rcRead += 'A';
		} else {
			rcRead += read[i];
		}
	}
	return rcRead;
}

char reverseComplementBase(const char& c) {
	if (c == 'A') {
		return 'T';
	} else if (c == 'C') {
		return 'G';
	} else if (c == 'G') {
		return 'C';
	} else if (c == 'T') {
		return 'A';
	} else {
		return c;
	}
}

void ReadWithAlignments::extractErrors(const seqan::Dna5String &genome) {
	// ignore searching errors in unmapped reads
	if (hasFlagUnmapped(records[0])) {
		return;
	}

	// ignore searching errors in soft-clipped non-chimeric reads
	/*if (records.size() == 1) {
	 for (size_t i = 0; i < length(records[0].cigar); ++i) {
	 if (records[0].cigar[i].operation == 'S') {
	 return;
	 }
	 }
	 }

	 // ignore searching errors in reads with chimeric breaks
	 if (records.size() > 1) {
	 throw std::runtime_error("records.size() > 1");
	 return;
	 }*/

	std::vector<CorrectionAligned> corrections;

	for (seqan::BamAlignmentRecord record : records) {

		std::string cigarString = "";

		// Extract CIGAR string and total read length (including unmapped regions)
		seqan::String< seqan::CigarElement<char> > cigar = record.cigar;
		unsigned readLength = 0;

		for (unsigned i = 0; i < length(cigar); ++i) {
			if (cigar[i].operation != 'D') {
				readLength += cigar[i].count;
			}
			cigarString += cigar[i].operation;
			cigarString += ":";
			cigarString += std::to_string(cigar[i].count);
			cigarString += " ";
		}

		assert(readLength == length(alignedRead.read));

		unsigned positionInRead = 0;
		unsigned hardClippedBases = 0;
		unsigned softClippedBases = 0;
		unsigned insertedBases = 0;
		unsigned deletedBases = 0;

		unsigned genomeBeginPos = record.beginPos;

		for (unsigned i = 0; i < length(cigar); ++i) {
			if (cigar[i].operation == 'H') { // hard clipping
				unsigned nucleotidePositionInReference = record.beginPos + positionInRead - insertedBases
						- softClippedBases + deletedBases;
				unsigned realPositionInRead = positionInRead + hardClippedBases; // the position of the first base of the chimeric break

				assert(realPositionInRead < readLength);

				if (realPositionInRead == 0) { // in this case, take the position of the last base of the chimeric break
					realPositionInRead += cigar[i].count;
				}

				// TODO: What if multiple chimeric breaks happen in a read? ... It should be fine, too. As long as there aren't multiple hard clipped regions in a single record.

				hardClippedBases += cigar[i].count;
				if (hasFlagRC(record)) {
					realPositionInRead = readLength - realPositionInRead - 1;
				}

				correctedRead.applyCorrection(
						CorrectionAligned(nucleotidePositionInReference,
								Correction(realPositionInRead, correctedRead.originalPositions[realPositionInRead], "",
										"", ErrorType::MULTIDEL))); // because we treat chimeric breaks as deletion of multiple bases
			} else if (cigar[i].operation == 'S') { // soft clipping
				for (size_t j = 0; j < cigar[i].count; ++j) {
					correctedRead.correctedRead.sequence[positionInRead + hardClippedBases + j] = 'S';
				}
				positionInRead += cigar[i].count;
				softClippedBases += cigar[i].count;

				//genomeBeginPos -= cigar[i].count;

				//throw std::runtime_error("Encountered soft clipping... this is currently a problem since we need the ground truth");

			} else if (cigar[i].operation == 'M') { // match ... check for substitutions later
				positionInRead += cigar[i].count;

			} else if (cigar[i].operation == 'I') { // insertion
				for (size_t j = 0; j < cigar[i].count; ++j) {
					unsigned nucleotidePositionRead = positionInRead;
					char nucleotideInRead = correctedRead.correctedRead.sequence[nucleotidePositionRead];

					if (nucleotideInRead == 'S') {
						throw std::runtime_error("this should not happen");
					}

					unsigned nucleotidePositionInReference = record.beginPos + nucleotidePositionRead - insertedBases
							- softClippedBases + deletedBases;
					size_t realPositionInRead = positionInRead + hardClippedBases;

					if (hasFlagRC(record)) {
						//realPositionInRead = readLength - realPositionInRead - 1;
						nucleotideInRead = reverseComplementBase(nucleotideInRead);
					}

					assert(realPositionInRead < readLength);

					std::string fromBase = "";
					fromBase += nucleotideInRead;

					size_t correctionPosition = realPositionInRead;
					if (hasFlagRC(record)) {
						correctionPosition = correctedRead.correctedRead.sequence.size() - correctionPosition - 1;
					}

					correctedRead.applyCorrection(
							CorrectionAligned(nucleotidePositionInReference,
									Correction(correctionPosition, correctedRead.originalPositions[correctionPosition],
											fromBase, "", ErrorType::INSERTION)));

				}
				insertedBases += cigar[i].count;

			} else if (cigar[i].operation == 'D') { // deletion
				size_t realPositionInRead = positionInRead + hardClippedBases - 1;

				/*if (hasFlagRC(record)) {
				 std::cout << "REVERSE COMPLEMENTED READ\n";

				 std::string genomeSubstr = "";
				 for (size_t i = record.beginPos; i < record.beginPos + readLength; ++i) {
				 genomeSubstr += genome[i];
				 }

				 std::cout << "genome substring:\n";
				 std::cout << genomeSubstr << "\n";

				 std::cout << "reverse complemented read before fixing deletion:\n";
				 std::cout << reverseComplementString(correctedRead.correctedRead.sequence) << "\n";
				 }

				 std::cout << "before fixing deletion:\n";
				 std::cout << correctedRead.correctedRead.sequence << "\n";*/

				assert(realPositionInRead < readLength);

				std::string fromBase = "";
				fromBase += correctedRead.correctedRead.sequence[realPositionInRead];

				if (fromBase == "S") {
					std::cout << cigarString << "\n";
					std::cout << "realpositionInRead:" << realPositionInRead << "\n";
					std::cout << correctedRead.correctedRead.sequence << "\n";
					throw std::runtime_error("this should not happen");
				}

				std::string toBases = "";

				size_t correctionPosition = realPositionInRead;
				if (hasFlagRC(record)) {
					correctionPosition = correctedRead.correctedRead.sequence.size() - correctionPosition - 1;
				}

				unsigned nucleotidePositionInReference = (positionInRead + hardClippedBases) + record.beginPos
						- softClippedBases;
				/*if (hasFlagRC(record)) {
				 nucleotidePositionInReference = record.beginPos - (positionInRead + hardClippedBases);
				 }*/

				if (nucleotidePositionInReference >= length(genome)) {
					std::string cigarString = "";
					for (unsigned i = 0; i < length(cigar); ++i) {
						cigarString += cigar[i].operation;
						cigarString += ":";
						cigarString += std::to_string(cigar[i].count);
						cigarString += " ";
					}
					std::cout << "cigar String:" << cigarString << "\n";
					std::cout << "beginPos: " << record.beginPos << "\n";
					std::cout << "read length: " << length(record.seq) << "\n";
					std::cout << "genome size: " << length(genome) << "\n";

					throw std::runtime_error("nucleotidePositionInReference >= length(genome)");
				}

				seqan::Dna5 nucleotideInReference;
				for (size_t j = 0; j < cigar[i].count; ++j) {
					deletedBases++;
					nucleotideInReference = genome[nucleotidePositionInReference + j];
					if (hasFlagRC(record)) {
						nucleotideInReference = reverseComplementBase(nucleotideInReference);
					}
					toBases += nucleotideInReference;
				}
				assert(toBases.size() == cigar[i].count);
				if (cigar[i].count == 1) {
					if (nucleotideInReference == 'A') {
						correctedRead.applyCorrection(
								CorrectionAligned(nucleotidePositionInReference,
										Correction(correctionPosition,
												correctedRead.originalPositions[correctionPosition], fromBase,
												fromBase + "A", ErrorType::DEL_OF_A)));
					} else if (nucleotideInReference == 'C') {
						correctedRead.applyCorrection(
								CorrectionAligned(nucleotidePositionInReference,
										Correction(correctionPosition,
												correctedRead.originalPositions[correctionPosition], fromBase,
												fromBase + "C", ErrorType::DEL_OF_C)));
					} else if (nucleotideInReference == 'G') {
						correctedRead.applyCorrection(
								CorrectionAligned(nucleotidePositionInReference,
										Correction(correctionPosition,
												correctedRead.originalPositions[correctionPosition], fromBase,
												fromBase + "G", ErrorType::DEL_OF_G)));
					} else if (nucleotideInReference == 'T') {
						correctedRead.applyCorrection(
								CorrectionAligned(nucleotidePositionInReference,
										Correction(correctionPosition,
												correctedRead.originalPositions[correctionPosition], fromBase,
												fromBase + "T", ErrorType::DEL_OF_T)));
					} else if (nucleotideInReference != 'N') {
						throw std::runtime_error("weird base in reference genome");
					}
				} else {
					correctedRead.applyCorrection(
							CorrectionAligned(nucleotidePositionInReference,
									Correction(correctionPosition, correctedRead.originalPositions[correctionPosition],
											fromBase, fromBase + toBases, ErrorType::MULTIDEL)));
				}
				positionInRead += cigar[i].count;

				/*std::cout << "after fixing deletion:\n";
				 std::cout << correctedRead.correctedRead.sequence << "\n";


				 std::cout << "reverse complemented read after fixing deletion:\n";
				 std::cout << reverseComplementString(correctedRead.correctedRead.sequence) << "\n";*/

			} else if (cigar[i].operation == 'P' || cigar[i].operation == 'N') {
				std::cout << "ERROR! There are letters I don't understand yet!" << cigar[i].operation << std::endl;
			}
		}

		// now that indels have been fixed, fix the substitution errors.

		/*std::cout << "before fixing substitution errors:\n";
		 std::cout << correctedRead.correctedRead.sequence << "\n";*/

		int genomeIdx = genomeBeginPos - 1;
		for (size_t i = 0; i < correctedRead.correctedRead.sequence.size(); ++i) {
			if (correctedRead.correctedRead.sequence[i] == 'S') {
				continue;
			} else {
				genomeIdx++;
			}

			std::string baseInRead = "";
			baseInRead += correctedRead.correctedRead.sequence[i];
			std::string baseInGenome = "";
			baseInGenome += genome[genomeIdx];

			size_t correctionPos = i;

			if (hasFlagRC(record)) {
				//std::cout << "REVERSE COMPLEMENTED READ\n";

				correctionPos = correctedRead.correctedRead.sequence.size() - i - 1;

				//std::string revComp = reverseComplementString(correctedRead.correctedRead.sequence);
				baseInRead = "";
				baseInRead += reverseComplementBase(correctedRead.correctedRead.sequence[correctionPos]);

				//baseInRead += revComp[i];
			}

			if (baseInRead != baseInGenome) {
				std::string debugString = "";
				for (size_t t = 0; t < readLength; ++t) {
					debugString += genome[genomeBeginPos + t];
				}

				ErrorType type;
				if (baseInGenome == "A") {
					type = ErrorType::SUB_FROM_A;

					if (hasFlagRC(record)) {
						type = ErrorType::SUB_FROM_T; // TODO: Check me.
					}

				} else if (baseInGenome == "C") {
					type = ErrorType::SUB_FROM_C;

					if (hasFlagRC(record)) {
						type = ErrorType::SUB_FROM_G; // TODO: Check me.
					}

				} else if (baseInGenome == "G") {
					type = ErrorType::SUB_FROM_G;

					if (hasFlagRC(record)) {
						type = ErrorType::SUB_FROM_C; // TODO: Check me.
					}

				} else if (baseInGenome == "T") {
					type = ErrorType::SUB_FROM_T;

					if (hasFlagRC(record)) {
						type = ErrorType::SUB_FROM_A; // TODO: Check me.
					}

				} else if (baseInGenome == "N") {
					throw std::runtime_error("base in genome is N");
					continue;
				} else {
					throw std::runtime_error("ERROR");
				}
				correctedRead.applyCorrection(
						CorrectionAligned(genomeIdx,
								Correction(correctionPos, correctedRead.originalPositions[correctionPos], baseInRead,
										baseInGenome, type)));
			}
		}

		/*std::cout << "after fixing substitution errors:\n";
		 std::cout << correctedRead.correctedRead.sequence << "\n";

		 std::cout << "beginPos: " << correctedRead.beginPos << "\n";

		 std::cout << cigarString << "\n";

		 MDTag mdTag;
		 extractMDTag(record, mdTag);
		 std::cout << mdTagToString(mdTag) << "\n";*/

		correctedRead.endPos = genomeIdx - 1;
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
