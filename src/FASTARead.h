/*
 * FASTARead.h
 *
 *  Created on: Jan 31, 2017
 *      Author: sarah
 */

#pragma once

#include <iostream>
#include <string>
#include "external/cereal/types/string.hpp"

class FASTARead {
public:
	FASTARead();
	FASTARead(const std::string &name, const std::string &sequence);

	std::string id;
	std::string seq;

	template<class Archive>
	void serialize(Archive & archive) {
		archive(id, seq); // serialize things by passing them to the archive
	}
};

inline std::ostream& operator<<(std::ostream & os, const FASTARead &read) {
	os << ">" << read.id << "\n" << read.seq;
	return os;
}

inline std::istream& operator>>(std::istream & is, FASTARead &read) {
	is >> read.id >> read.seq;
	read.id = read.id.substr(1, read.id.size() - 1);
	return is;
}
