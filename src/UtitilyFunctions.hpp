/*
 * UtitilyFunctions.hpp
 * This source file contains widely used utility functions.
 *
 *  Created on: Feb 6, 2017
 *      Author: sarah
 */

#pragma once

#include <cstddef>
#include <unordered_map>
#include <utility>

struct pairhash {
public:
  template <typename T, typename U>
  std::size_t operator()(const std::pair<T, U> &x) const
  {
    return 13 * std::hash<T>()(x.first) + 37 * std::hash<U>()(x.second);
  }
};

template <typename T, typename U>
inline std::ostream& operator<<(std::ostream& os, const std::pair<T, U> &x)
{
	return os << x.first << ' ' << x.second;
}

template <typename T, typename U>
inline std::istream& operator>>(std::istream& is, std::pair<T, U> &x)
{
	return is >> x.first >> x.second;
}
