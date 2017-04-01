/*
 * MotifTreeEntry.hpp
 *
 *  Created on: Feb 7, 2017
 *      Author: sarah
 */

#pragma once

#include <stddef.h>
#include <unordered_map>
#include <cereal/types/unordered_map.hpp>

#include "../../ErrorType.h"

class MotifTreeEntry {
public:
	MotifTreeEntry();
	void reset();

	std::unordered_map<ErrorType, size_t> numErrors;
	std::unordered_map<ErrorType, double> zScore;

	template<class Archive>
	void serialize(Archive & archive) {
		archive(numErrors, zScore); // serialize things by passing them to the archive
	}
};
