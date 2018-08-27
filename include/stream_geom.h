//! \file stream_geom.h
//! For describing STREAM geometries
#ifndef STREAM_STREAM_GEOM_H
#define STREAM_STREAM_GEOM_H

namespace stream {

//! Describes an (x,y,z) coordinate in 3D space.
struct Position {
  double x;  //!< x-coordinate
  double y;  //!< y-coordinate
  double z;  //!< z-coordinate
};

} // namespace stream

#endif // STREAM_STREAM_GEOM_H
