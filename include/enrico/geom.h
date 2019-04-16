//! \file geom.h
//! For describing ENRICO geometries
#ifndef ENRICO_GEOM_H
#define ENRICO_GEOM_H

namespace enrico {

//! Describes an (x,y,z) coordinate in 3D space.
struct Position {
  double x; //!< x-coordinate
  double y; //!< y-coordinate
  double z; //!< z-coordinate
};

} // namespace enrico

#endif // ENRICO_GEOM_H
