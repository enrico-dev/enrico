# =============================================================================
#
# mkuserfile.cmake: Generates a .f file from a .usr file.
# 
# A direct port of Nek5000/core/mkuserfile. This is intended to be the COMMAND 
# for a call to add_custom_command() that generates the .usr file. It can also 
# be run directly from the command-line:
#     $ cmake -DCASENAME=casename [-DCVODE] [-DCMT] -P mkuserfile.cmake
#
# =============================================================================

# =============================================================================
# Read .usr file
# =============================================================================

file(READ ${INFILE_DIR}/${CASENAME}.usr usr_str)
string(TOLOWER "${usr_str}" usr_str_lower)

# =============================================================================
# Fail if userq is defined.  It must be specifically defined for coupling
# =============================================================================

if(usr_str_lower MATCHES "subroutine.*userq")
  message(FATAL_ERROR "userq() was defined in ${CASENAME}.usr. For coupling, userq() is \
  predifined and you cannot redefine your own.  Remove the definition of userq() in \
  ${CASENAME}.usr and run make again.")
endif()

# =============================================================================
# Add standard subroutines
# =============================================================================

if(NOT usr_str_lower MATCHES "subroutine.*usrsetvert")
  string(CONCAT usr_str "${usr_str}" 
"
c automatically added by cmake
      subroutine usrsetvert(glo_num,nel,nx,ny,nz) ! to modify glo_num
      integer*8 glo_num(1)

      return
      end
"
  )
endif()

if(NOT usr_str_lower MATCHES "subroutine.*userqtl")
  string(CONCAT usr_str "${usr_str}" 
"
c automatically added by cmake
      subroutine userqtl

      call userqtl_scig

      return
      end
"
  )
endif()

# =============================================================================
# Add CVODE subroutines
# =============================================================================

if(CVODE)

  if(NOT usr_str_lower MATCHES "include.*cvode_aux.*h")
    string(CONCAT usr_str "${usr_str}" 
"
c automatically added by makenek
\#include \"cvode_aux.h\"
"
    )
  endif()

  if(NOT usr_str_lower MATCHES "include.*cvode_jtimes.*h")
    string(CONCAT usr_str "${usr_str}" 
"
c automatically added by makenek
\#include \"cvode_jtimes.h\"
"
    )
  endif()

  if(NOT usr_str_lower MATCHES "include.*cvode_preco.*h")
    string(CONCAT usr_str "${usr_str}" 
"
c automatically added by makenek
\#include \"cvode_preco_dummy.h\"
"
    )
  endif()

endif()

# =============================================================================
# Add CMT subroutines
# =============================================================================

if(CMT)

  if(NOT usr_str_lower MATCHES "subroutine.*cmt_usrflt")
    string(CONCAT usr_str "${usr_str}" 
"
c automatically added by makenek
      subroutine cmt_usrflt(rmult) ! user defined filter
      include 'SIZE'
      real rmult(lx1)
      call rone(rmult,lx1)
      return
      end
"
    )
  endif()

  if(NOT usr_str_lower MATCHES "subroutine.*cmt_userflux")
    string(CONCAT usr_str "${usr_str}" 
"
c automatically added by makenek
      subroutine cmt_userflux ! user defined flux
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      include 'CMTDATA'
      real fluxout(lx1*lz1)
      return
      end
"
    )
  endif()

  if(NOT usr_str_lower MATCHES "subroutine.*cmt_usereos")
    string(CONCAT usr_str "${usr_str}" 
"
c automatically added by makenek
      subroutine cmt_userEOS ! user defined EOS 
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      include 'CMTDATA'

      return
      end
"
    )
  endif()

  # TODO: Define a full path relative to cmake build dir
  if(EXISTS "cmtparticles.usrp")
    message(STATUS "Particles found CMT")
  else()
    string(CONCAT usr_str "${usr_str}" 
"
c automatically added by makenek
      subroutine usr_particles_init ! used for particles
      return
      end
c
c automatically added by makenek
      subroutine usr_particles_solver ! used for particles
      return
      end
c
c automatically added by makenek
      subroutine usr_particles_io(istep) ! used for particles
      integer istep
      return
      end
"
    )
  endif()

endif() 

# =============================================================================
# Write .f file
# =============================================================================

file(WRITE ${OUTFILE_DIR}/${CASENAME}.f ${usr_str})
