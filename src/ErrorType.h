/*
 * ErrorType.h
 *
 *  Created on: Jan 31, 2017
 *      Author: sarah
 */

#pragma once

#include <type_traits>
#include <iostream>

enum class ErrorType : unsigned {
	CORRECT = 0, INSERTION = 1, SUB_FROM_A = 2, SUB_FROM_C = 3, SUB_FROM_G = 4, SUB_FROM_T = 5, DEL_OF_A = 6, DEL_OF_C = 7, DEL_OF_G = 8, DEL_OF_T = 9, MULTIDEL = 10, NODEL = 11
};

inline std::string errorTypeToString(const ErrorType &type) {
	std::string res = "";
	switch (type) {
	case ErrorType::CORRECT:
		res = "CORRECT";
		break;
	case ErrorType::INSERTION:
		res = "INSERTION";
		break;
	case ErrorType::SUB_FROM_A:
		res = "SUB_FROM_A";
		break;
	case ErrorType::SUB_FROM_C:
		res = "SUB_FROM_C";
		break;
	case ErrorType::SUB_FROM_G:
		res = "SUB_FROM_G";
		break;
	case ErrorType::SUB_FROM_T:
		res = "SUB_FROM_T";
		break;
	case ErrorType::DEL_OF_A:
		res = "DEL_OF_A";
		break;
	case ErrorType::DEL_OF_C:
		res = "DEL_OF_C";
		break;
	case ErrorType::DEL_OF_G:
		res = "DEL_OF_G";
		break;
	case ErrorType::DEL_OF_T:
		res = "DEL_OF_T";
		break;
	case ErrorType::MULTIDEL:
		res = "MULTIDEL";
		break;
	case ErrorType::NODEL:
		res = "NODEL";
		break;
	};
	return res;
}

inline unsigned errorTypeToNumber(const ErrorType &type) {
	return (unsigned) type;
}

inline ErrorType errorTypeFromString(const std::string &typeString) {
	if (typeString == "CORRECT") {
		return ErrorType::CORRECT;
	} else if (typeString == "INSERTION") {
		return ErrorType::INSERTION;
	} else if (typeString == "SUB_FROM_A") {
		return ErrorType::SUB_FROM_A;
	} else if (typeString == "SUB_FROM_C") {
		return ErrorType::SUB_FROM_C;
	} else if (typeString == "SUB_FROM_G") {
		return ErrorType::SUB_FROM_G;
	} else if (typeString == "SUB_FROM_T") {
		return ErrorType::SUB_FROM_T;
	} else if (typeString == "DEL_OF_A") {
		return ErrorType::DEL_OF_A;
	} else if (typeString == "DEL_OF_C") {
		return ErrorType::DEL_OF_C;
	} else if (typeString == "DEL_OF_G") {
		return ErrorType::DEL_OF_G;
	} else if (typeString == "DEL_OF_T") {
		return ErrorType::DEL_OF_T;
	} else if (typeString == "MULTIDEL") {
		return ErrorType::MULTIDEL;
	} else if (typeString == "NODEL") {
		return ErrorType::NODEL;
	}
	throw std::runtime_error(typeString + " does not describe a valid error type!");
}

inline std::ostream& operator<<(std::ostream & os, const ErrorType &type) {
	return os << errorTypeToString(type);
}

inline std::istream& operator>>(std::istream &is, ErrorType &type) {
	std::string errorTypeString;
	is >> errorTypeString;
	type = errorTypeFromString(errorTypeString);
	return is;
}

// code from http://stackoverflow.com/questions/261963/how-can-i-iterate-over-an-enum
template<typename C, C beginVal, C endVal>
class Iterator {
	typedef typename std::underlying_type<C>::type val_t;
	int val;
public:
	Iterator(const C & f) :
			val(static_cast<val_t>(f)) {
	}
	Iterator() :
			val(static_cast<val_t>(beginVal)) {
	}
	Iterator operator++() {
		++val;
		return *this;
	}
	C operator*() {
		return static_cast<C>(val);
	}
	Iterator begin() {
		return *this;
	} //default ctor is good
	Iterator end() {
		static const Iterator endIter = ++Iterator(endVal); // cache it
		return endIter;
	}
	bool operator!=(const Iterator& i) {
		return val != i.val;
	}
};

namespace cereal
{
  template <class Archive> inline
  std::string save_minimal( Archive const &, ErrorType const & t )
  {
    return errorTypeToString( t );
  }

  template <class Archive> inline
  void load_minimal( Archive const &, ErrorType & t, std::string const & value )
  {
    t = errorTypeFromString( value );
  }
}

typedef Iterator<ErrorType, ErrorType::CORRECT, ErrorType::NODEL> errorTypeIterator;
typedef Iterator<ErrorType, ErrorType::CORRECT, ErrorType::SUB_FROM_T> errorTypesCurrentBase;
typedef Iterator<ErrorType, ErrorType::DEL_OF_A, ErrorType::NODEL> errorTypesNextGap;
typedef Iterator<ErrorType, ErrorType::INSERTION, ErrorType::MULTIDEL> errorTypesError;
