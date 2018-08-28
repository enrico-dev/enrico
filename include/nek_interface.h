#ifndef STREAM_NEK_INTERFACE_H
#define STREAM_NEK_INTERFACE_H

#include "nek_mangling.h"
#include "stream_geom.h"

extern "C" {

// From libnek5000
void C2F_nek_init(const int* intracomm);
void C2F_nek_end();
void C2F_nek_solve();

int nek_get_global_elem_centroid(int global_elem, stream::Position* centroid);
int nek_get_global_elem(int local_elem);
int nek_get_local_elem(int global_elem);
int nek_get_lelg();
int nek_get_lelt();
int nek_get_lx1();
int nek_get_nelgt();
int nek_get_nelt();

void nek_init_step();
void nek_step();
void nek_finalize_step();
};

#endif // STREAM_NEK_INTERFACE_H
