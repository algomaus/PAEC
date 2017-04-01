/*
 * FASTQRead.cpp
 *
 *  Created on: Jan 31, 2017
 *      Author: sarah
 */

#include "FASTQRead.h"

#include <cassert>
#include <iostream>

FASTQRead::FASTQRead() {
	id = "";
	sequence = "";
	quality = "";
}

FASTQRead::FASTQRead(const std::string &name, const std::string &readSequence, const std::string &qual) {
	id = name;
	sequence = readSequence;
	quality = qual;
}

FASTQRead::FASTQRead(const FASTARead &fasta, const std::string &qual) {
	id = fasta.id;
	sequence = fasta.seq;
	assert(qual.size() == sequence.size());
	quality = qual;
}

