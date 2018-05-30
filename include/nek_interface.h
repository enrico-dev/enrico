#ifndef STREAM_NEK_INTERFACE_H
#define STREAM_NEK_INTERFACE_H

#include "nek_mangling.h"

extern "C" {

// From libnek5000
void C2F_nek_init(const int *intracomm);
void C2F_nek_init_step();
void C2F_nek_step();
void C2F_nek_finalize_step();
void C2F_nek_end();
void C2F_nek_solve();

// From nek_interface.f90
int nek_get_local_el_centroid(int local_el, double *x, double *y, double *z);

};

#endif //STREAM_NEK_INTERFACE_H
