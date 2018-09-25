#ifndef STREAM_ERROR_H
#define STREAM_ERROR_H
#include <stdexcept>
#include <string>
#include "openmc/capi.h"

namespace stream {

constexpr int E_SUCCESS = 0;

//! Wrapper for error checking a called function
//!
//! This follows OpenMC's conventions that an error is < 0 (as opposed to != 0).
//!
//! \param err The error code returned by the called function.
//! \param msg The message passed displayed if an error occurs.
inline void err_chk(const int err, const char* msg)
{
  if (err < E_SUCCESS)
    throw std::runtime_error(msg);
}

inline void err_chk(const int err, const std::string msg)
{
  if (err < E_SUCCESS)
    throw std::runtime_error(msg);
}

inline void err_chk(const int err)
{
  if (err < E_SUCCESS)
    throw std::runtime_error(openmc_err_msg);
}

}
#endif //STREAM_ERROR_H

