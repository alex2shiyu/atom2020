!=========================================================================!
! project : rambutan
! program : context module
! history : Apr 26, 2011
! authors : xidai and duliang {duleung@gmail.com}
! purpose : important arrays defined for main program 
! comment :
!=========================================================================!

!-------------------------------------------------------------------------!
!>>> basis states related matrix
!-------------------------------------------------------------------------!
  module atomic_basis

     implicit none

! number of dimensions of each (isub) subspace
     integer, public, allocatable, save :: nstat(:)

! decimal representation of basis state in Fock space
     integer, public, allocatable, save :: basis(:)

! serial number of a decimal representated configuration
     integer, public, allocatable, save :: invsn(:)

! binary representation of a decimal represented basis state
     integer, public, allocatable, save :: invcd(:, :)

  end module atomic_basis

!-------------------------------------------------------------------------!
!>>> angular momentum related matrix
!-------------------------------------------------------------------------!
  module atomic_ammat
     use constants

     implicit none
! added by sypeng
! sy operator matrix in corresponding basis(see details below)
     complex(dp), public, allocatable, save :: c4mat_i(:, :)
     complex(dp), public, allocatable, save :: c4mat_r(:, :)
     complex(dp), public, allocatable, save :: c4mat_t(:, :)
     complex(dp), public, allocatable, save :: c4mat_n(:, :)
     complex(dp), public, allocatable, save :: c4mat_bd(:, :)
     complex(dp), public, allocatable, save :: c4mat_irrep(:, :)
     complex(dp), public, allocatable, save :: c2mat_r(:, :)
     complex(dp), public, allocatable, save :: c2mat_t(:, :)
     complex(dp), public, allocatable, save :: c2mat_n(:, :)
     complex(dp), public, allocatable, save :: c2mat_bd(:, :)
     complex(dp), public, allocatable, save :: c2mat_irrep(:, :)
     complex(dp), public, allocatable, save :: sgmdmat_r(:, :)
     complex(dp), public, allocatable, save :: sgmdmat_t(:, :)
     complex(dp), public, allocatable, save :: sgmdmat_n(:, :)
     complex(dp), public, allocatable, save :: sgmdmat_bd(:, :)
     complex(dp), public, allocatable, save :: sgmdmat_irrep(:, :)
     complex(dp), public, allocatable, save :: sgmvmat_r(:, :)
     complex(dp), public, allocatable, save :: sgmvmat_t(:, :)
     complex(dp), public, allocatable, save :: sgmvmat_n(:, :)
     complex(dp), public, allocatable, save :: sgmvmat_bd(:, :)
     complex(dp), public, allocatable, save :: sgmvmat_irrep(:, :)
! eigenvalue of c4 operator in many body basis
     real(dp),    public, allocatable, save :: c4eig(:)
     integer,     public, allocatable, save :: c4degs(:)
     integer,     public, allocatable, save :: c4virrepdegs(:)
     complex(dp), public, allocatable, save :: c4v_umt(:, :)
     complex(dp), public, allocatable, save :: c4v_rdct_umt(:, :)
     integer,     public, allocatable, save :: c4v_umt_1d(:)
     integer,     public, allocatable, save :: irrep_flag(:)
     integer,     public, allocatable, save :: irrep_type(:)
     integer,     public, allocatable, save :: num_irrep(:)
     integer,     public,              save :: num_rdu
     integer,     public,              save :: num_2d
! ----------------------------------------------------------
! ss matrix in many body basis
     complex(dp), public, allocatable, save :: ssmat_n(:, :)

! ll matrix in many body basis
     complex(dp), public, allocatable, save :: llmat_n(:, :)

! jj matrix in many body basis
     complex(dp), public, allocatable, save :: jjmat_n(:, :)

! lx operator matrix in {lz, sz} basis
     complex(dp), public, allocatable, save :: lxmat_i(:, :)
! ly operator matrix in {lz, sz} basis
     complex(dp), public, allocatable, save :: lymat_i(:, :)
! lz operator matrix in {lz, sz} basis
     complex(dp), public, allocatable, save :: lzmat_i(:, :)
! sx operator matrix in {lz, sz} basis
     complex(dp), public, allocatable, save :: sxmat_i(:, :)
! sy operator matrix in {lz, sz} basis
     complex(dp), public, allocatable, save :: symat_i(:, :)
! sz operator matrix in {lz, sz} basis
     complex(dp), public, allocatable, save :: szmat_i(:, :)
! jx operator matrix in {lz, sz} basis
     complex(dp), public, allocatable, save :: jxmat_i(:, :)
! jy operator matrix in {lz, sz} basis
     complex(dp), public, allocatable, save :: jymat_i(:, :)
! jz operator matrix in {lz, sz} basis
     complex(dp), public, allocatable, save :: jzmat_i(:, :)

! lx operator matrix in {j2, jz} basis
     complex(dp), public, allocatable, save :: lxmat_r(:, :)
! lm operator matrix in {j2, jz} basis
     complex(dp), public, allocatable, save :: lymat_r(:, :)
! lz operator matrix in {j2, jz} basis
     complex(dp), public, allocatable, save :: lzmat_r(:, :)
! sp operator matrix in {j2, jz} basis
     complex(dp), public, allocatable, save :: sxmat_r(:, :)
! sm operator matrix in {j2, jz} basis
     complex(dp), public, allocatable, save :: symat_r(:, :)
! sz operator matrix in {j2, jz} basis
     complex(dp), public, allocatable, save :: szmat_r(:, :)
! jp operator matrix in {j2, jz} basis
     complex(dp), public, allocatable, save :: jxmat_r(:, :)
! jm operator matrix in {j2, jz} basis
     complex(dp), public, allocatable, save :: jymat_r(:, :)
! jz operator matrix in {j2, jz} basis
     complex(dp), public, allocatable, save :: jzmat_r(:, :)

! lp operator matrix in {j2, jz} basis
     complex(dp), public, allocatable, save :: lxmat_t(:, :)
! lm operator matrix in {j2, jz} basis
     complex(dp), public, allocatable, save :: lymat_t(:, :)
! lz operator matrix in {j2, jz} basis
     complex(dp), public, allocatable, save :: lzmat_t(:, :)
! sp operator matrix in {j2, jz} basis
     complex(dp), public, allocatable, save :: sxmat_t(:, :)
! sm operator matrix in {j2, jz} basis
     complex(dp), public, allocatable, save :: symat_t(:, :)
! sz operator matrix in {j2, jz} basis
     complex(dp), public, allocatable, save :: szmat_t(:, :)
! jp operator matrix in {j2, jz} basis
     complex(dp), public, allocatable, save :: jxmat_t(:, :)
! jm operator matrix in {j2, jz} basis
     complex(dp), public, allocatable, save :: jymat_t(:, :)
! jz operator matrix in {j2, jz} basis
     complex(dp), public, allocatable, save :: jzmat_t(:, :)

! lp operator matrix in many body basis
     complex(dp), public, allocatable, save :: lxmat_n(:, :)
! lm operator matrix in many body basis
     complex(dp), public, allocatable, save :: lymat_n(:, :)
! lz operator matrix in many body basis
     complex(dp), public, allocatable, save :: lzmat_n(:, :)
! sp operator matrix in many body basis
     complex(dp), public, allocatable, save :: sxmat_n(:, :)
! sm operator matrix in many body basis
     complex(dp), public, allocatable, save :: symat_n(:, :)
! sz operator matrix in many body basis
     complex(dp), public, allocatable, save :: szmat_n(:, :)
! jp operator matrix in many body basis
     complex(dp), public, allocatable, save :: jxmat_n(:, :)
! jm operator matrix in many body basis
     complex(dp), public, allocatable, save :: jymat_n(:, :)
! jz operator matrix in many body basis
     complex(dp), public, allocatable, save :: jzmat_n(:, :)


! sz operator matrix in many body basis
     complex(dp), public, allocatable, save :: szmat_irrep(:, :)
     complex(dp), public, allocatable, save :: ssmat_irrep(:, :)
     complex(dp), public, allocatable, save :: ssmat_bd(:, :)




! eigenvalue of j2 operator in many body basis
     real(dp), public, allocatable, save :: jjeig(:)

! eigenvalue of jz operator in many body basis
     real(dp), public, allocatable, save :: jzeig(:)

! eigenvalue of s2 operator in many body basis
     real(dp), public, allocatable, save :: sseig(:)

! eigenvalue of sz operator in many body basis
     real(dp), public, allocatable, save :: szeig(:)

! eigenvalue of c4 operator in mang body basis
!     complex(dp), public, allocatable, save :: c4eig(:)
  end module atomic_ammat

!-------------------------------------------------------------------------!
!>>> Hamiltonian martix related variables
!-------------------------------------------------------------------------!
  module atomic_hmtrx
     use constants

     implicit none

! gutzwiller variational parameters matrix
     integer, public, allocatable, save :: Hdegs(:)
     integer, public, allocatable, save :: vpmmat(:, :)

! eigenvalues of atomic Hamiltonian matrix
     real(dp), public, allocatable, save :: heigs(:)
     real(dp), public, allocatable, save :: elocs(:)
     real(dp), public, allocatable, save :: neigs(:)

! eigenvector of atomic Hamiltonian matrix
     complex(dp), public, allocatable, save :: heigv(:, :)
    
! impurity energy level
     complex(dp), public, allocatable, save :: cemat(:, :)

! atomic Hamiltonian matrix
     complex(dp), public, allocatable, save :: hmat(:, :)

! transformation matrix from {px, sz} basis to {lz, sz} basis
     complex(dp), public, allocatable, save :: dmtrx(:, :)

! transformation matrix from {lz, sz} basis to {j2, jz} basis
     complex(dp), public, allocatable, save :: cgmat(:, :)

! transformation matrix from {px, sz} basis to {j2, jz} basis
     complex(dp), public, allocatable, save :: amtrx(:, :)

! spin orbit coupling matrix in {lz, sz} single particle basis
     complex(dp), public, allocatable, save :: somat(:, :)

! spin orbit coupling matrix in {j2, jz} single particle basis
     complex(dp), public, allocatable, save :: somat_t(:, :)

! auxiliary complex(dp) matrix for temperary use
     complex(dp), public, allocatable, save :: zauxs(:, :)

! matrix form of construct operator in many body basis
     complex(dp), public, allocatable, save :: fpadd(:, :, :)

! matrix form of eliminate operator in many body basis
     complex(dp), public, allocatable, save :: fpcut(:, :, :)

! coefficents matrix for generalized interaction U in {px,sz} basis
     complex(dp), public, allocatable, save :: cumat(:, :, :, :)

! coefficents matrix for generalized interaction U in {j2, jz} basis
     complex(dp), public, allocatable, save :: cumat_t(:, :, :, :)

! nn operator matrix in many body basis
     complex(dp), public, allocatable, save :: npmat(:, :)

! \sum_\alpha(n_{\alpha\minus} - n_{\alpha\plus} )^2 operator matrix in many body basis
! which is actually single occupied number
     complex(dp), public, allocatable, save :: sgnpmat(:, :)
     complex(dp), public, allocatable, save :: sgnpmat_irrep(:, :)
! nn operator matrix in many body basis
     complex(dp), public, allocatable, save :: npmat_irrep(:, :)
! atomic Hamiltonian matrix in c4v's irrep basis
     complex(dp), public, allocatable, save :: hmat_irrep(:, :)


! eigenvalue of n_{sg} operator in many body basis
     real(dp), public, allocatable, save :: sgeigs(:)

  end module atomic_hmtrx

!-------------------------------------------------------------------------!
!>>> memory managment subroutines (allocate and deallocate memory)
!-------------------------------------------------------------------------!
  module context
     use constants
     use control

     use atomic_ammat
     use atomic_basis
     use atomic_hmtrx

     implicit none

! status flag
     integer, private :: istat

! declaration of module procedures: allocate memory
     public :: atomic_allocate_memory_basis
     public :: atomic_allocate_memory_hmtrx
     public :: atomic_allocate_memory_ammat

! declaration of module procedures: deallocate memory
     public :: atomic_deallocate_memory_basis
     public :: atomic_deallocate_memory_hmtrx

     contains

!>>> allocate memory atomic_basis
     subroutine atomic_allocate_memory_basis()

         implicit none

         allocate(nstat(0:norbs), stat=istat)
         allocate(basis(1:ncfgs), stat=istat)

         allocate(invsn(0:2**norbs-1), stat=istat)

         allocate(invcd(1:norbs, 1:ncfgs), stat=istat)

         if (istat /= 0) then
             stop "error happened in atomic_allocate_memory_basis"
         endif ! back if (istat /= 0) block

! initialize the variables
         nstat = 0; basis = 0
         invsn = 0; invcd = 0

         return
     end subroutine atomic_allocate_memory_basis

!>>> allocate memory atomic_angle
     subroutine atomic_allocate_memory_ammat()
         implicit none
! added by sypeng
         allocate(c4mat_r(norbs, norbs), stat=istat)
         allocate(c4mat_i(norbs, norbs), stat=istat)
         allocate(c4mat_t(norbs, norbs), stat=istat)
         allocate(c4mat_n(ncfgs, ncfgs), stat=istat)
         allocate(c4mat_bd(ncfgs, ncfgs), stat=istat)
         allocate(c4mat_irrep(ncfgs, ncfgs), stat=istat)
         allocate(hmat_irrep(ncfgs, ncfgs), stat=istat)
         c4mat_r = czero;     c4mat_i = czero;   c4mat_t = czero
         c4mat_n = czero;     c4mat_bd = czero
         c4mat_irrep = czero; hmat_irrep = czero
         allocate(c4eig(ncfgs), stat=istat); c4eig = czero
!        --------------         
         allocate(c2mat_r(norbs, norbs), stat=istat)
         allocate(c2mat_t(norbs, norbs), stat=istat)
         allocate(c2mat_n(ncfgs, ncfgs), stat=istat)
         allocate(c2mat_bd(ncfgs, ncfgs), stat=istat)
         allocate(c2mat_irrep(ncfgs, ncfgs), stat=istat)
         c2mat_r = czero;  c2mat_t = czero; c2mat_n = czero
         c2mat_bd = czero; c2mat_irrep = czero
         allocate(sgmvmat_r(norbs, norbs), stat=istat)
         allocate(sgmvmat_t(norbs, norbs), stat=istat)
         allocate(sgmvmat_n(ncfgs, ncfgs), stat=istat)
         allocate(sgmvmat_bd(ncfgs, ncfgs), stat=istat)
         allocate(sgmvmat_irrep(ncfgs, ncfgs), stat=istat)
         sgmvmat_r = czero; sgmvmat_t = czero; sgmvmat_n = czero
         sgmvmat_bd = czero; sgmvmat_irrep = czero
         allocate(sgmdmat_r(norbs, norbs), stat=istat)
         allocate(sgmdmat_t(norbs, norbs), stat=istat)
         allocate(sgmdmat_n(ncfgs, ncfgs), stat=istat)
         allocate(sgmdmat_bd(ncfgs, ncfgs), stat=istat)
         allocate(sgmdmat_irrep(ncfgs, ncfgs), stat=istat)
         sgmdmat_r = czero; sgmdmat_t = czero; sgmdmat_n = czero
         sgmdmat_bd = czero; sgmdmat_irrep = czero
         allocate(c4degs(ncfgs), stat=istat)
         allocate(c4virrepdegs(ncfgs), stat=istat)
         allocate(irrep_flag(ncfgs), stat=istat)
         allocate(irrep_type(ncfgs), stat=istat)
         allocate(num_irrep(6), stat=istat)
         allocate(c4v_umt(ncfgs, ncfgs), stat=istat)
         allocate(c4v_rdct_umt(2, 2), stat=istat)
         allocate(c4v_umt_1d(ncfgs), stat=istat)
         allocate(sgeigs(ncfgs), stat=istat)
         c4eig      = 0.d0;  c4degs       = 0; c4v_umt      = dcmplx(0.d0, 0.d0)
         irrep_flag = 0   ;  irrep_type   = 0; c4v_umt_1d   = 0
         num_rdu    = 0   ;  num_2d       = 0; c4v_rdct_umt = dcmplx(0.d0, 0.d0)
         num_irrep  = 0   ;  c4virrepdegs = 0; sgeigs       = 0.d0
! --------------------------------------------------

         allocate(sxmat_r(norbs, norbs), stat=istat)
         allocate(sxmat_i(norbs, norbs), stat=istat)
         allocate(sxmat_t(norbs, norbs), stat=istat)
         allocate(sxmat_n(ncfgs, ncfgs), stat=istat)
         sxmat_r = czero; sxmat_i = czero; sxmat_t = czero; sxmat_n = czero

         allocate(symat_r(norbs, norbs), stat=istat)
         allocate(symat_i(norbs, norbs), stat=istat)
         allocate(symat_t(norbs, norbs), stat=istat)
         allocate(symat_n(ncfgs, ncfgs), stat=istat)
         symat_r = czero; symat_i = czero; symat_t = czero; symat_n = czero

         allocate(szmat_r(norbs, norbs), stat=istat)
         allocate(szmat_i(norbs, norbs), stat=istat)
         allocate(szmat_t(norbs, norbs), stat=istat)
         allocate(szmat_n(ncfgs, ncfgs), stat=istat)
         szmat_r = czero; szmat_i = czero; szmat_t = czero; szmat_n = czero

         allocate(lxmat_r(norbs, norbs), stat=istat)
         allocate(lxmat_i(norbs, norbs), stat=istat)
         allocate(lxmat_t(norbs, norbs), stat=istat)
         allocate(lxmat_n(ncfgs, ncfgs), stat=istat)
         lxmat_r = czero; lxmat_i = czero; lxmat_t = czero; lxmat_n = czero

         allocate(lymat_r(norbs, norbs), stat=istat)
         allocate(lymat_i(norbs, norbs), stat=istat)
         allocate(lymat_t(norbs, norbs), stat=istat)
         allocate(lymat_n(ncfgs, ncfgs), stat=istat)
         lymat_r = czero; lymat_i = czero; lymat_t = czero; lymat_n = czero

         allocate(lzmat_r(norbs, norbs), stat=istat)
         allocate(lzmat_i(norbs, norbs), stat=istat)
         allocate(lzmat_t(norbs, norbs), stat=istat)
         allocate(lzmat_n(ncfgs, ncfgs), stat=istat)
         lzmat_r = czero; lzmat_i = czero; lzmat_t = czero; lzmat_n = czero

         allocate(jxmat_r(norbs, norbs), stat=istat)
         allocate(jxmat_i(norbs, norbs), stat=istat)
         allocate(jxmat_t(norbs, norbs), stat=istat)
         allocate(jxmat_n(ncfgs, ncfgs), stat=istat)
         jxmat_r = czero; jxmat_i = czero; jxmat_t = czero; jxmat_n = czero

         allocate(jymat_r(norbs, norbs), stat=istat)
         allocate(jymat_i(norbs, norbs), stat=istat)
         allocate(jymat_t(norbs, norbs), stat=istat)
         allocate(jymat_n(ncfgs, ncfgs), stat=istat)
         jymat_r = czero; jymat_i = czero; jymat_t = czero; jymat_n = czero

         allocate(jzmat_r(norbs, norbs), stat=istat)
         allocate(jzmat_i(norbs, norbs), stat=istat)
         allocate(jzmat_t(norbs, norbs), stat=istat)
         allocate(jzmat_n(ncfgs, ncfgs), stat=istat)
         jzmat_r = czero; jzmat_i = czero; jzmat_t = czero; jzmat_n = czero

         allocate(jjeig(ncfgs), stat=istat); jjeig = dzero
         allocate(jzeig(ncfgs), stat=istat); jzeig = dzero
         allocate(sseig(ncfgs), stat=istat); sseig = dzero
         allocate(szeig(ncfgs), stat=istat); szeig = dzero
!         allocate(c4eig(ncfgs), stat=istat); c4eig = czero
         allocate(ssmat_n(ncfgs, ncfgs), stat=istat); ssmat_n = czero
         allocate(llmat_n(ncfgs, ncfgs), stat=istat); llmat_n = czero
         allocate(jjmat_n(ncfgs, ncfgs), stat=istat); jjmat_n = czero

         allocate(ssmat_irrep(ncfgs, ncfgs), stat=istat); ssmat_irrep = czero
         allocate(szmat_irrep(ncfgs, ncfgs), stat=istat); szmat_irrep = czero
         allocate(ssmat_bd(ncfgs, ncfgs), stat=istat); ssmat_bd = czero
         return
     end subroutine atomic_allocate_memory_ammat

!>>> allocate memory atomic_hmtrx
     subroutine atomic_allocate_memory_hmtrx()
         implicit none

         allocate(vpmmat(ncfgs, ncfgs), stat=istat)

         allocate(elocs(norbs), stat=istat)
         allocate(neigs(ncfgs), stat=istat)
         allocate(heigs(ncfgs), stat=istat)
         allocate(hdegs(ncfgs), stat=istat)
         allocate(heigv(ncfgs, ncfgs), stat=istat)
         allocate(zauxs(ncfgs, ncfgs), stat=istat)

         allocate(hmat(ncfgs, ncfgs), stat=istat)
         allocate(cemat(norbs, norbs), stat=istat)

         allocate(  dmtrx(norbs, norbs), stat=istat)
         allocate(  cgmat(norbs, norbs), stat=istat)
         allocate(  amtrx(norbs, norbs), stat=istat)
         allocate(  somat(norbs, norbs), stat=istat)
         allocate(somat_t(norbs, norbs), stat=istat)

         allocate(  npmat(ncfgs, ncfgs), stat=istat)
         allocate(npmat_irrep(ncfgs, ncfgs), stat=istat)
         allocate(sgnpmat(ncfgs, ncfgs), stat=istat)
         allocate(sgnpmat_irrep(ncfgs, ncfgs), stat=istat)

         allocate(fpadd(ncfgs, ncfgs, norbs), stat=istat)
         allocate(fpcut(ncfgs, ncfgs, norbs), stat=istat)

         allocate(  cumat(norbs, norbs, norbs, norbs), stat=istat)
         allocate(cumat_t(norbs, norbs, norbs, norbs), stat=istat)

         if (istat /= 0) then
             stop "error happened in atomic_allocate_memory_hmtrx"
         endif ! back if (istat /= 0) block

! initialize the variables
         vpmmat = 0

         cemat  = czero; heigs   = dzero
         hmat  = czero; heigv   = czero

         fpadd = czero; fpcut = czero

         cgmat = czero
         dmtrx = czero; amtrx   = czero
         cumat = czero; cumat_t = czero
         somat = czero; somat_t = czero

         return
     end subroutine atomic_allocate_memory_hmtrx

!>>> deallocate memory atomic_basis
     subroutine atomic_deallocate_memory_basis()
         implicit none

         if (allocated(nstat)) deallocate(nstat)
         if (allocated(invsn)) deallocate(invsn)
         if (allocated(invcd)) deallocate(invcd)
         if (allocated(basis)) deallocate(basis)

         return
     end subroutine atomic_deallocate_memory_basis

!>>> deallocate memory atomic_hmtrx
     subroutine atomic_deallocate_memory_hmtrx()
         implicit none

         if (allocated(cemat)) deallocate(cemat)
         if (allocated(hmat)) deallocate(hmat)

         if (allocated(cumat)) deallocate(cumat)
         if (allocated(somat)) deallocate(somat)

         if (allocated(heigs)) deallocate(heigs)
         if (allocated(heigv)) deallocate(heigv)

         if (allocated(cumat_t)) deallocate(cumat_t)
         if (allocated(somat_t)) deallocate(somat_t)
         
         return
     end subroutine atomic_deallocate_memory_hmtrx

  end module context
