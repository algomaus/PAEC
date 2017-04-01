/*
 * MotifTreeNode.hpp
 *
 *  Created on: Feb 7, 2017
 *      Author: sarah
 */

#pragma once

#include <stddef.h>
#include <vector>
#include <cereal/types/vector.hpp>

#include "MotifTreeEntry.h"

class MotifTreeNode {
public:
	MotifTreeNode();
	void reset();

	char base;
	std::vector<MotifTreeEntry> entries; // size of the vector = level in the tree = length of the motif

	template<class Archive>
	void serialize(Archive & archive) {
		archive(base, entries); // serialize things by passing them to the archive
	}
};


