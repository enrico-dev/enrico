#ifndef STREAM_NEK_INTERFACE_H
#define STREAM_NEK_INTERFACE_H

#include "nek_mangling.h"
#include "stream_geom.h"

extern "C" {

// From libnek5000
void C2F_nek_init(const int *intracomm);
void C2F_nek_end();
void C2F_nek_solve();

/** Retrieves an array of centriods for a given array of local element numbers.
 *
 * This is the interoperable declaration for the Fortran function \ref nek_get_lelt_centroids
 *
 * @param[in]  lelts   An array of local element numbers.
 * @param[in]  n_lelts The number of local elements in *lelts*.
 * @param[out] ctroids An array of centroids corresponding to each local element in *lelts*.
 */
inline int nek_get_lelt_centroids(const int *lelts, const int n_lelts,
                                  Position *ctroids);

void nek_init_step();
void nek_step();
void nek_finalize_step();
};

#endif // STREAM_NEK_INTERFACE_H
