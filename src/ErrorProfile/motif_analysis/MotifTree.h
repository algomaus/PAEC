/*
 * MotifTree.hpp
 *
 *  Created on: Feb 7, 2017
 *      Author: sarah
 */

#pragma once

#include <unordered_map>
#include <cmath>
#include <vector>
#include <string>
#include "../../ErrorType.h"

#include "MotifTreeNode.h"

class MotifTree {
public:
	MotifTree();
	MotifTree(const int &maxMotifSize);
	MotifTreeNode& operator[](std::size_t idx)       { return nodes[idx]; }
	const MotifTreeNode& operator[](std::size_t idx) const { return nodes[idx]; }
	MotifTreeNode& operator[](const std::string &motifString)       { return nodes[motifToIndex(motifString)]; }
	const MotifTreeNode& operator[](const std::string &motifString) const { return nodes[motifToIndex(motifString)]; }

	size_t motifToIndex(const std::string &motif) const;
	std::string indexToMotif(size_t index);
	size_t size();
	void reset();

	template<class Archive>
	void serialize(Archive & archive) {
		archive(_maxMotifSize, nodes); // serialize things by passing them to the archive
	}
private:
	int _maxMotifSize;
	std::vector<MotifTreeNode> nodes;
};
