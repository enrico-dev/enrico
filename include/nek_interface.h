#ifndef C2F_HEADER_INCLUDED
#define C2F_HEADER_INCLUDED

/* Mangling for Fortran global symbols without underscores. */
#define C2F_GLOBAL(name,NAME) name##_

/* Mangling for Fortran global symbols with underscores. */
#define C2F_GLOBAL_(name,NAME) name##_

/* Mangling for Fortran module symbols without underscores. */
#define C2F_MODULE(mod_name,name, mod_NAME,NAME) __##mod_name##_MOD_##name

/* Mangling for Fortran module symbols with underscores. */
#define C2F_MODULE_(mod_name,name, mod_NAME,NAME) __##mod_name##_MOD_##name

/*--------------------------------------------------------------------------*/
/* Mangle some symbols automatically.                                       */
#define C2F_nek_init C2F_GLOBAL_(nek_init, NEK_INIT)
#define C2F_nek_init_step C2F_GLOBAL_(nek_init_step, NEK_INIT_STEP)
#define C2F_nek_step C2F_GLOBAL_(nek_step, NEK_STEP)
#define C2F_nek_finalize_step C2F_GLOBAL_(nek_finalize_step, NEK_FINALIZE_STEP)
#define C2F_nek_end C2F_GLOBAL_(nek_end, NEK_END)

#endif
