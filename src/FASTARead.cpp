/*
 * FASTARead.cpp
 *
 *  Created on: Jan 31, 2017
 *      Author: sarah
 */

#include "FASTARead.h"

FASTARead::FASTARead() {
	id = "";
	seq = "";
}

FASTARead::FASTARead(const std::string &name, const std::string &sequence) {
	id = name;
	seq = sequence;
}
