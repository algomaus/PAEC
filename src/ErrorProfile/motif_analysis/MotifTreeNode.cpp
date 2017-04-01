/*
 * MotifTreeNode.cpp
 *
 *  Created on: Feb 7, 2017
 *      Author: sarah
 */

#include "MotifTreeNode.h"

MotifTreeNode::MotifTreeNode() {
	base = '_';
}

void MotifTreeNode::reset() {
	for (size_t i = 0; i < entries.size(); ++i) {
		entries[i].reset();
	}
}
