#ifndef STREAM_NEK_INTERFACE_H
#define STREAM_NEK_INTERFACE_H

#include "nek_mangling.h"
#include "stream_geom.h"

extern "C" {

// From libnek5000
void C2F_nek_init(const int *intracomm);
void C2F_nek_init_step();
void C2F_nek_step();
void C2F_nek_finalize_step();
void C2F_nek_end();
void C2F_nek_solve();

// From nek_interface.f90
inline int nek_get_lelt_centroids(const int *lelts, const int n_lelts,
                                  Position *ctroids);
};

#endif // STREAM_NEK_INTERFACE_H
