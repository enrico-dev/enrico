#ifndef ENRICO_UTILS_H
#define ENRICO_UTILS_H

#include <cmath>

namespace enrico {

//! Compare two numbers for
inline bool soft_equiv(double value, double reference)
{
  constexpr double rel_precision = 1e-12;
  constexpr double abs_threshold = 1e-14;

  // Typical case: relative error comparison to reference
  if (std::abs(value - reference) < rel_precision * std::abs(reference)) {
    return true;
  }

  // If one is within the absolute threshold of zero, and the other within
  // relative of zero, they're equal
  if ((std::abs(reference) < abs_threshold) && (std::abs(value) < rel_precision)) {
    return true;
  }
  return (std::abs(value) < abs_threshold) && (std::abs(reference) < rel_precision);
}

} // namespace enrico

#endif // ENRICO_UTILS_H
