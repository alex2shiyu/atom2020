!-------------------------------------------------------------------------!
! project : rambutan
! program : atomic_print_header
!           atomic_print_footer
!           atomic_print_summary
! author  : duliang (email:duleung@gamil.com)
! history : 11/28/2011
! purpose : print heading, ending, and summary of current project
! comment :
!-------------------------------------------------------------------------!
!>>> print the heading information
  subroutine atomic_print_header()
     use constants
     use control

     implicit none

     write(mystd, '(2x,a)') ''
     write(mystd, '(2x,a)') '>>> rotational invariant gutzwiller method'
     write(mystd, *)

     write(mystd, '(2x,a)') 'version: '
     write(mystd, '(2x,a)') 'develop: '
     write(mystd, '(2x,a)') 'support: '
     write(mystd, *)

     write(mystd, '(2x,a)') 'rambutan >>> running'

# if defined (MPI)

     write(mystd, '(2x,a,i3)') 'rambutan >>> parallel: Y >>> nprocs:', nprocs

# else   /* MPI */

     write(mystd, '(2x,a,i3)') 'rambutan >>> parallel: N >>> nprocs:', 1

# endif  /* MPI */

     write(mystd, *)

     return
  end subroutine atomic_print_header

!>>> print the ending information
  subroutine atomic_print_footer()
     use constants
     use control

     implicit none

! used to record the time information
     real(dp) :: time

! obtain time information
     call cpu_time(time)

     write(mystd, '(2X,A,F10.2,A)') 'rambutan >>> total time spent:', time, 's'
     write(mystd, *)

     write(mystd, '(2X,A)') 'rambutan >>> hope you good luck. Bye!'
     write(mystd, '(2X,A)') 'rambutan >>> ending'

     return
  end subroutine atomic_print_footer

!>>> print the running parameters
  subroutine atomic_print_summary()
     use constants
     use control

     implicit none

     write(mystd, '(2X, A)') 'pasture >>> parameters list:'

     write(mystd, '(2(4X, A, I10))')   'nband :', nband , 'nspin :', nspin
     write(mystd, '(2(4X, A, I10))')   'norbs :', norbs , 'ncfgs :', ncfgs

     write(mystd, '(2(4X, A, F10.5))') 'Uc    :', Uc    , 'Uv    :', Uv
     write(mystd, '(2(4X, A, F10.5))') 'Jz    :', Jz    , 'Js    :', Js
     write(mystd, '(2(4X, A, F10.5))') 'Jp    :', Jp    , 'mune  :', mune

     write(mystd, *)

     return
  end subroutine atomic_print_summary
