/*
 * KmerType.h
 *
 *  Created on: Jan 31, 2017
 *      Author: sarah
 */

#pragma once

enum class KmerType {
	REPEAT = 0, TRUSTED = 1, UNTRUSTED = 2
};

inline std::string kmerTypeToString(const KmerType &type) {
	std::string res = "";
	switch (type) {
	case KmerType::REPEAT:
		res = "REPEAT";
		break;
	case KmerType::TRUSTED:
		res = "TRUSTED";
		break;
	case KmerType::UNTRUSTED:
		res = "UNTRUSTED";
		break;
	};
	return res;
}

inline KmerType kmerTypeFromString(const std::string &typeString) {
	if (typeString == "REPEAT") {
		return KmerType::REPEAT;
	} else if (typeString == "TRUSTED") {
		return KmerType::TRUSTED;
	} else if (typeString == "UNTRUSTED") {
		return KmerType::UNTRUSTED;
	} else {
		throw std::runtime_error(typeString + " is not a valid KmerType!");
	}
}

inline unsigned kmerTypeToNumber(const KmerType &type) {
	unsigned res = 0;
	switch (type) {
	case KmerType::REPEAT:
		res = 0;
		break;
	case KmerType::TRUSTED:
		res = 1;
		break;
	case KmerType::UNTRUSTED:
		res = 2;
		break;
	};
	return res;
}

inline KmerType kmerTypeFromNumber(const int &num) {
	KmerType res;
	switch (num) {
	case 0:
		res = KmerType::REPEAT;
		break;
	case 1:
		res = KmerType::TRUSTED;
		break;
	case 2:
		res = KmerType::UNTRUSTED;
		break;
	default:
		throw std::runtime_error("Invalid KmerType number: " + std::to_string(num));
	};
	return res;
}

inline std::ostream& operator<<(std::ostream & os, KmerType type) {
	return os << kmerTypeToString(type);
}

namespace cereal {
template<class Archive> inline std::string save_minimal(Archive const &, KmerType const & t) {
	return kmerTypeToString(t);
}

template<class Archive> inline
void load_minimal(Archive const &, KmerType & t, std::string const & value) {
	t = kmerTypeFromString(value);
}
}
