/*
 * FASTQRead.h
 *
 *  Created on: Jan 31, 2017
 *      Author: sarah
 */

#pragma once

#include <iostream>
#include <string>
#include "external/cereal/types/string.hpp"

#include "FASTARead.h"

class FASTQRead {
public:
	FASTQRead();
	FASTQRead(const std::string &name, const std::string &readSequence, const std::string &qual);
	FASTQRead(const FASTARead &fasta, const std::string &qual);

	std::string id;
	std::string sequence;
	std::string quality;

	template<class Archive>
	void serialize(Archive & archive) {
		archive(id, sequence, quality); // serialize things by passing them to the archive
	}
};

inline std::ostream& operator<<(std::ostream & os, const FASTQRead &read) {
	os << read.id << "\n" << read.sequence << "\n+\n" << read.quality;
	return os;
}

inline std::istream& operator>>(std::istream & is, FASTQRead &read) {
	std::string garbage;
	is >> read.id >> read.sequence >> garbage >> read.quality;
	read.id = read.id.substr(1, read.id.size() - 1);
	return is;
}
