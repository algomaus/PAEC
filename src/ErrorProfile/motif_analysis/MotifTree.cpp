/*
 * MotifTree.cpp
 *
 *  Created on: Feb 7, 2017
 *      Author: sarah
 */

#include "MotifTree.h"

MotifTree::MotifTree() {
	_maxMotifSize = 0;
}
MotifTree::MotifTree(const int &maxMotifSize) {
	_maxMotifSize = maxMotifSize;
	size_t numMotifs = 0;
	for (int l = 1; l <= _maxMotifSize; ++l) {
		numMotifs += std::pow(5, l);
	}
	//std::cout << "Creating motif tree for " << numMotifs << " motifs...\n";
	nodes.resize(numMotifs);
	for (int l = 1; l <= _maxMotifSize; ++l) {
		size_t lengthOffset = 0;
		for (int j = 1; j < l; ++j) {
			lengthOffset += std::pow(5, j);
		}

		for (size_t i = 0; i < std::pow(5, l); ++i) {
			nodes[lengthOffset + i].entries.resize(l);
			switch (i % 5) {
			case 0:
				nodes[lengthOffset + i].base = 'A';
				break;
			case 1:
				nodes[lengthOffset + i].base = 'C';
				break;
			case 2:
				nodes[lengthOffset + i].base = 'G';
				break;
			case 3:
				nodes[lengthOffset + i].base = 'T';
				break;
			case 4:
				nodes[lengthOffset + i].base = 'N';
				break;
			}
		}
	}
	//std::cout << "Finished creation of motif tree. \n";
}

size_t MotifTree::motifToIndex(const std::string &motif) const {
	size_t index = 0;
	size_t l = motif.size();
	size_t lengthOffset = 0;
	for (size_t j = 1; j < l; ++j) {
		lengthOffset += std::pow(5, j);
	}
	int i = l - 1;
	while (i >= 0) {
		size_t base_i = std::pow(5, l - i - 1);
		if (motif[i] == 'A') {
			index += 0 * base_i;
		} else if (motif[i] == 'C') {
			index += 1 * base_i;
		} else if (motif[i] == 'G') {
			index += 2 * base_i;
		} else if (motif[i] == 'T') {
			index += 3 * base_i;
		} else if (motif[i] == 'N') {
			index += 4 * base_i;
		}
		i--;
	}
	return lengthOffset + index;
}

std::string MotifTree::indexToMotif(size_t index) {
	size_t l;
	size_t old = 0;
	if (index < 5) {
		l = 1;
	} else {
		l = 2;
		old = 5;
		while (index >= old + std::pow(5, l)) {
			old += std::pow(5, l);
			l++;
		}
	}

	// now l is the size of the motif.
	size_t lengthOffset = old;
	index = index - lengthOffset;

	std::string motifString = "";
	while (true) {
		size_t remainder = index % 5;
		index /= 5;
		switch (remainder) {
		case 0:
			motifString = "A" + motifString;
			break;
		case 1:
			motifString = "C" + motifString;
			break;
		case 2:
			motifString = "G" + motifString;
			break;
		case 3:
			motifString = "T" + motifString;
			break;
		case 4:
			motifString = "N" + motifString;
			break;
		}
		if (index == 0)
			break;
	}

	// fill the remaining ones with "A"...
	while (motifString.size() < l) {
		motifString = "A" + motifString;
	}

	return motifString;
}

size_t MotifTree::size() {
	return nodes.size();
}

void MotifTree::reset() {
	for (size_t i = 0; i < nodes.size(); ++i) {
		nodes[i].reset();
	}
}
