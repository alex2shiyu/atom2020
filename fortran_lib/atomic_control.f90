!=========================================================================!
! project : rambutan
! program : control
! history : Apr 26, 2011
! author  : duliang (duleung@gmail.com)
! purpose : some impartant common control variables
! comment : 
!=========================================================================!
  module control
     use constants

     implicit none

!-------------------------------------------------------------------------!
!>>> atomic Hamiltonian parameters
!-------------------------------------------------------------------------!
! number of bands
     integer, public, save :: nband

! number of spins
     integer, public, save :: nspin

! number of orbits 
     integer, public, save :: norbs

! number of configurations
     integer, public, save :: ncfgs

! number of vpms
     integer, public, save :: nvpms

! number of total electrons
     integer, public, save :: ntots

! intraorbital Coulomb interaction
     real(dp), public, save :: Uc  

! interorbital Coulomb interaction
     real(dp), public, save :: Uv  

! Hund's exchange interaction
     real(dp), public, save :: Jz  

! spin-flip interaction
     real(dp), public, save :: Js

! pair-hopping interaction
     real(dp), public, save :: Jp

! spin-orbit coupling interaction
     real(dp), public, save :: lamb

! hopping integral between nearest neighbor sites
     real(dp), public, save :: thop

! chemical potential
     real(dp), public, save :: mune

! inverse temperture (pseudo temperture)
     real(dp), public, save :: beta

!-------------------------------------------------------------------------!
!>>> MPI related common variables
!-------------------------------------------------------------------------!
! number of processors
     integer, public, save :: nprocs

! the rank of the controller process
     integer, public, save :: master

! the rank of the current process
     integer, public, save :: myrank

  end module control
