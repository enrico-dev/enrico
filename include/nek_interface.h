#ifndef STREAM_NEK_INTERFACE_H
#define STREAM_NEK_INTERFACE_H

#include "nek_mangling.h"

extern "C" {

void C2F_nek_init(const int *intracomm);
void C2F_nek_init_step();
void C2F_nek_step();
void C2F_nek_finalize_step();
void C2F_nek_end();
void C2F_nek_solve();

};

#endif //STREAM_NEK_INTERFACE_H
