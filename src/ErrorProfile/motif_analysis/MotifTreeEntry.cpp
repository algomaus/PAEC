/*
 * MotifTreeEntry.cpp
 *
 *  Created on: Feb 7, 2017
 *      Author: sarah
 */

#include "MotifTreeEntry.h"

MotifTreeEntry::MotifTreeEntry() {
	for (ErrorType type : errorTypesError()) {
		numErrors[type] = 0;
		zScore[type] = 0;
	}
}

void MotifTreeEntry::reset() {
	for (ErrorType type : errorTypesError()) {
		numErrors[type] = 0;
		zScore[type] = 0;
	}
}
