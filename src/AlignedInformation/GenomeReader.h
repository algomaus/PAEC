/*
 * GenomeReader.h
 *
 *  Created on: Jan 31, 2017
 *      Author: sarah
 */

#pragma once

#include <seqan/sequence.h>
#include <string>

class GenomeReader {
public:
	GenomeReader();
	static std::string readGenome(const std::string &genomeFilename);
};
