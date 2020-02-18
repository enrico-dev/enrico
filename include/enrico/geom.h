//! \file geom.h
//! For describing ENRICO geometries
#ifndef ENRICO_GEOM_H
#define ENRICO_GEOM_H

namespace enrico {

//! Describes an (x,y,z) coordinate in 3D space.
struct Position {
  // Constructors
  Position() {}
  Position(double x0, double y0, double z0)
    : x{x0}
    , y{y0}
    , z{z0}
  {}

  double x{}; //!< x-coordinate
  double y{}; //!< y-coordinate
  double z{}; //!< z-coordinate
};

} // namespace enrico

#endif // ENRICO_GEOM_H
