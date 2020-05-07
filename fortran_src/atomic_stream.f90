!-------------------------------------------------------------------------!
! project : rambutan
! program : atomic_config
! history : 01/02/2012
! author  : duliang (duleung@gmail.com)
! purpose : setup atomic Hamiltonian parameters
! comment : 
!-------------------------------------------------------------------------!
!>>> atomic hamiltonian parameters
  subroutine atomic_config()
     use constants
     use control

     implicit none

! local variables
! check whether the input file (dft.atom.in) exist
     logical :: exists

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

     nband = 3            ! number of bands
     nspin = 2            ! number of spins
     norbs = 6            ! number of orbits
     ncfgs = 20           ! number of configurations

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

     Uc    = 4.00_dp      ! intraorbital Coulomb interaction
     Uv    = 3.00_dp      ! interorbital Coulomb interaction
     Jz    = 0.50_dp      ! Hund's exchange interaction
     Js    = 0.50_dp      ! spin-flip interaction
     Jp    = 0.50_dp      ! pair-hopping interaction
     lamb  = 0.00_dp      ! spin-orbit coupling parameter
     thop  = 0.50_dp      ! hopping integral parameter
     mune  = 0.00_dp      ! chemical potential
     beta  = 10.0_dp      ! inverse temperture (pseudo temperture)

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

! read from input file if it exists
     exists = .false.

! inquire the input file status
     inquire( file="dft.atom.in", exist=exists )

! read parameters from dft.atom.in
     if ( exists .eqv. .true. ) then
         open( mytmp, file="dft.atom.in", status="old" )

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
         read(mytmp, *)
         read(mytmp, *)
         read(mytmp, *)

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
         read(mytmp, *) nband
         read(mytmp, *) nspin
         read(mytmp, *) norbs
         read(mytmp, *) ntots
         read(mytmp, *) ncfgs

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
         read(mytmp, *)

         read(mytmp, *) Uc
         read(mytmp, *) Uv
         read(mytmp, *) Jz
         read(mytmp, *) Js
         read(mytmp, *) Jp

         read(mytmp, *) lamb
         read(mytmp, *) thop
         read(mytmp, *) mune
         read(mytmp, *) beta

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

         close(mytmp)
     endif ! back if ( exists .eqv. .true. ) block

     return
  end subroutine atomic_config

!>>> allocate memory and initialize them
  subroutine atomic_setup_array()
     use constants
     use control
     use context

     implicit none

! allocate memory for angular related matrix
     call atomic_allocate_memory_ammat()

! allocate memory for basis sheet
     call atomic_allocate_memory_basis()

! allocate memory for hamiltonian matrix
     call atomic_allocate_memory_hmtrx()

     return
  end subroutine atomic_setup_array

!>>> deallocate memory and finalize them
  subroutine atomic_final_array()
     use constants
     use control
     use context

     implicit none

! deallocate memory for basis sheet
     call atomic_deallocate_memory_basis()

! deallocate memory for hamiltonian matrix
     call atomic_deallocate_memory_hmtrx()

     return
  end subroutine atomic_final_array

!>>> calculate combination algebra 
  function state_pick(ntiny, nlarg) result(value)
     implicit none

! external variables
     integer, intent(in) :: ntiny
     integer, intent(in) :: nlarg

! local variables
     integer :: i

! auxiliary integer variable
     integer :: nlow

! numberator of the combination algebra
     real(8) :: numer

! denominator of the combination algebra
     real(8) :: denom

! result value of the combination algebra
     integer :: value

! transform the combination algebra
     nlow = min(ntiny, nlarg-ntiny)

! numerator in combination algebra
     numer = 1.d0
     do i=nlarg-nlow+1,nlarg
        numer = numer * dble(i)
     enddo ! over i={nlarg-nlow+1,nlarg} loop

! denominator in combination algebra
     denom = 1.d0
     do i=1,nlow
        denom = denom * dble(i)
     enddo ! over i={1,nlow} loop

! result value
     value = nint(numer / denom)

     return
  end function state_pick

!-------------------------------------------------------------------------!
!>>> build key matrix in atomic problem
!-------------------------------------------------------------------------!
  subroutine atomic_jovian(nstat, dmtrx, cgmat, cemat, somat, cumat)
     use constants
     use control

     implicit none

! external arguments
! external functions
     integer, external :: state_pick

! number of configuration in each subspace
     integer, intent(out) :: nstat(0:norbs)

! impurity energy level
     complex(dp), intent(out) :: cemat(norbs, norbs)

! transformation matrix bwtween {px, sz} basis and {lz, sz} basis
     complex(dp), intent(out) :: dmtrx(norbs, norbs)

! transformation matrix between {lz, sz} and {j2, jz} basis
     complex(dp), intent(out) :: cgmat(norbs, norbs)

! spin orbit coupling matrix in {lz, sz} basis
     complex(dp), intent(out) :: somat(norbs, norbs)

! general interaction coefficents in {px, sz} basis
     complex(dp), intent(out) :: cumat(norbs, norbs, norbs, norbs)

! local variables
! loop index over good quantum number
     integer :: ibit

! loop index over orbits
     integer :: iorb
     integer :: jorb

! initialize dimension of each subspace
     do ibit=0,norbs
         nstat(ibit) = state_pick( ibit, norbs )
         print*, ibit, nstat(ibit)
     enddo ! over ibit={0,norbs} loop 

! build impurity energy level matrix
     call atomic_make_cemat(norbs, cemat)

! build transformation matrix from {lz, sz} to {j2, jz} basis
!     noted by syp, we don't consider soc here
     call atomic_make_cgmat(nband, nspin, norbs, cgmat)

! build spin-orbit coupling matrix in {lz, sz} basis
!
     call atomic_make_somatp(nband, nspin, norbs, lamb, somat)

! build transformation matrix from {px,sz} to {lz, sz} basis
     call atomic_make_dmatd(nband, nspin, norbs, dmtrx)

! build general interaction coefficents matrix in {lz, sz} basis
     call atomic_make_cumat( norbs, Uc, Uv, Jz, Js, Jp, cumat )

     return
  end subroutine atomic_jovian

!>>> impurity energy level
  subroutine atomic_make_cemat( norbs, cemat )
     use constants

     implicit none

! external arguments
! number of orbits
     integer, intent(in) :: norbs

! impurity energy level
     complex(dp), intent(out) :: cemat(norbs, norbs)

! local variable
! check whether the input file (rtgw.eimp.in) exist
     logical :: exists

! status while reading data
     integer :: ierr

! orbital index
     integer :: iorb
     integer :: jorb

! auxiliary real(dp) variables
     real(dp) :: xreal, ximag

! initialize eimp to be czero
     cemat = czero

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     do iorb=1,norbs
         if ( iorb .le. norbs-2 ) then
             cemat(iorb, iorb) = - 0.0d0
         else
             cemat(iorb, iorb) = + 0.0d0
         endif
     enddo ! over iorb={1,norbs} loop

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
! read from  input file if it exists
     exists = .false.

! inquire the input file status
     inquire( file="dft.eimp.in", exist=exists )

! read from rtgw.eimp.in
     if ( exists .eqv. .true. ) then
         open( mytmp, file="dft.eimp.in", status="old" )

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
         read(mytmp, *)
         read(mytmp, *)
         read(mytmp, *)

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
         do while ( .true. )
             read(mytmp, *, iostat=ierr) iorb, jorb, xreal, ximag
             if ( ierr /= 0 ) exit
             cemat(iorb, jorb) = dcmplx(xreal, ximag)
         enddo ! over while ( .true. ) loop

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

         close(mytmp)
         print *,"read dft.eimp.in successfully!"
     else
         cemat(1, 1)   = 0.133069
         cemat(2, 2)   = 0.133069
         cemat(3, 3)   = 0.133069
         cemat(4, 4)   = 0.133069
         cemat(5, 5)   = 0.298848
         cemat(6, 6)   = 0.298848
         cemat(7, 7)   = 0.267877
         cemat(8, 8)   = 0.267877
         cemat(9, 9)   = 0.004367
         cemat(10, 10) = 0.004367
         print *,"dft.eimp.in doesn't exist, be careful!"
     endif ! back if ( exists .eqv. .true. ) block

# if defined (Duliang)
     open(mytst, file='test-eimp.dat', status='unknown')
     do iorb=1,norbs
         do jorb=1,norbs
             if ( abs(cemat(iorb, jorb)) .lt. eps9 ) cycle
             write(mytst, '(2i8, 2f17.10)') iorb, jorb, cemat(iorb, jorb)
         enddo ! over jorb={1,norbs} loop
     enddo ! over iorb={1,norbs} loop
     close(mytst)
# endif /* Duliang */
     return
  end subroutine atomic_make_cemat

!>>>  $U_{\alpha\betta\delta\gamma}$ coefficent matrix for Uij  <<<!
!--------------------three band for example------------------------!
!> $f_{\alpha}^{\dagger}f_{\betta}^{\dagger}f_{\delta}f_{\gamma}$ <!
! norbs   bandindex(?band)    spinindex(?spin)    representation   !
!   1          1               1->\uparrow             1\up        !
!   2          1               0->\doarrow             1\do        !
!   3          2               1->\uparrow             2\up        !
!   4          2               0->\doarrow             2\do        !
!   5          3               1->\uparrow             3\up        !
!   6          3               0->\doarrow             3\do        !
!------------------------------------------------------------------!
  subroutine atomic_make_cumat( norbs, Uc, Uv, Jz, Js, Jp, cumat )
     use constants

     implicit none

! external arguments
! number of orbits
     integer, intent(in) :: norbs

! general coulomb interaction parameters
     real(dp), intent(in) :: Uc
     real(dp), intent(in) :: Uv
     real(dp), intent(in) :: Jz
     real(dp), intent(in) :: Js
     real(dp), intent(in) :: Jp

! general coulomb interaction matrix
     complex(dp), intent(out) :: cumat(norbs, norbs, norbs, norbs)

! local varibales
! loop index over orbits
     integer :: alpha, betta
     integer :: delta, gamma

! band index and spin index
! band index of alpha and betta
     integer :: aband, bband 

! band index of delta and gamma
     integer :: dband, gband 

! spin index of alpha and betta
     integer :: aspin, bspin 

! spin index of delta and gamma
     integer :: dspin, gspin 

! auxiliary variables
     real(dp) :: res

! initialize cumat to zero
     cumat = czero

! loop for creation operator
     alphaloop: do alpha=1,norbs-1
     bettaloop: do betta=alpha+1,norbs

! loop for destroy operator
        gammaloop: do gamma=1,norbs-1
        deltaloop: do delta=gamma+1,norbs
            aband = (alpha+1)/2; aspin = mod(alpha,2)
            bband = (betta+1)/2; bspin = mod(betta,2)
            gband = (gamma+1)/2; gspin = mod(gamma,2)
            dband = (delta+1)/2; dspin = mod(delta,2)

! here we use "res" due to overlap between "Uv and Jz"
            res = dzero

! intraorbital Coulomb interaction
            if ((alpha.eq.gamma) .and. (betta.eq.delta)) then
                if ((aband.eq.bband) .and. (aspin.ne.bspin)) then
                    res = res + Uc
!#                  write(mystd, "(a3,4i4)") " Uc", alpha, betta, delta, gamma
                endif
            endif

! interorbital Coulomb interaction
            if ((alpha.eq.gamma) .and. (betta.eq.delta)) then
                if (aband .ne. bband) then
                    res = res + Uv
!#                  write(mystd, "(a3,4i4)") " Uv", alpha, betta, delta, gamma
                endif
            endif

! Hund's exchange interaction 
            if ((alpha.eq.gamma) .and. (betta.eq.delta)) then
                if ((aband.ne.bband) .and. (aspin.eq.bspin)) then
                    res = res - Jz
!#                  write(mystd, "(a3,4i4)") " Jz", alpha, betta, delta, gamma
                endif
            endif
           
! spin flip term
            if ((aband.eq.gband) .and. (bband.eq.dband)) then
                if ((aspin.ne.gspin) .and. (bspin.ne.dspin) .and. (aspin.ne.bspin)) then
                    res = res - Js
!#                  write(mystd, "(a3,4i4)") " Js", alpha, betta, delta, gamma
                endif
            endif
          
! pair hopping term
            if ((aband.eq.bband) .and. (dband.eq.gband) .and. (aband.ne.dband)) then
                if ((aspin.ne.bspin) .and. (dspin.ne.gspin) .and. (aspin.eq.gspin)) then
                    res = res + Jp
!#                  write(mystd, "(a3,4i4)") " Jp", alpha, betta, delta, gamma
                endif
            endif
                 
            cumat(alpha, betta, delta, gamma) = res

        enddo deltaloop ! over delta={gamma+1,norbs} loop
        enddo gammaloop ! over gamma={1,norbs-1} loop
     enddo bettaloop ! over betta={alpha+1,norbs} loop
     enddo alphaloop ! over alpha={1,norbs-1} loop

# if defined (Duliang)
     open(mytst, file='test-cumat.dat', status='unknown')
     do alpha=1,norbs
     do betta=1,norbs
         do delta=1,norbs
         do gamma=1,norbs
             if ( abs(cumat(alpha, betta, delta, gamma)) .lt. eps6 ) cycle
             write(mytst, '(4i8, 2f17.10)') alpha, betta, delta, gamma, &
             cumat(alpha, betta, delta, gamma)
         enddo ! over gamma={1,norbs} loop
         enddo ! over delta={1,norbs} loop
     enddo ! over betta={1,norbs} loop
     enddo ! over alpha={1,norbs} loop
     close(mytst)
# endif /* Duliang */
     return
  end subroutine atomic_make_cumat

  subroutine atomic_tran_cumat(norbs, amtrx, cumat, cumat_t)
     use constants

     implicit none

! number of orbits
     integer, intent(in) :: norbs

! transformation matrix from {px, sz} basis to {j2, jz} basis
     complex(dp), intent(in) :: amtrx(norbs, norbs)

! coefficents matrix for generalized interaction U in {lz, sz} basis
     complex(dp), intent(in) :: cumat(norbs, norbs, norbs, norbs)
 
! coefficents matrix for generalized interaction U in {j2, jz} basis
     complex(dp), intent(out) :: cumat_t(norbs, norbs, norbs, norbs)

! local varoables
! loop index over orbits in {lz, sz} single particle basis
     integer :: alpha1, alpha2
     integer :: alpha3, alpha4

! loop index over orbits in {j2, jz} single particle basis
     integer :: sigma1, sigma2
     integer :: sigma3, sigma4

! auxiliary complex(dp) variables
     complex(dp) :: ctmp
     complex(dp) :: ztmp

! initialize cumat_t to be zero
     cumat_t = czero

     sigma1loop: do sigma1=1,norbs
     sigma2loop: do sigma2=1,norbs
     sigma3loop: do sigma3=1,norbs
     sigma4loop: do sigma4=1,norbs
         ctmp = czero

         alpha1loop: do alpha1=1,norbs
         alpha2loop: do alpha2=1,norbs
         alpha3loop: do alpha3=1,norbs
         alpha4loop: do alpha4=1,norbs
             if (abs(cumat(alpha1, alpha2, alpha3, alpha4)) .lt. eps6) cycle
             ctmp = ctmp + cumat(alpha1, alpha2, alpha3, alpha4)          &
                  * conjg(amtrx(alpha1, sigma1)) * amtrx(alpha3, sigma3)  &
                  * conjg(amtrx(alpha2, sigma2)) * amtrx(alpha4, sigma4)
         enddo alpha4loop ! over alpha4={1,norbs} loop
         enddo alpha3loop ! over alpha3={1,norbs} loop
         enddo alpha2loop ! over alpha2={1,norbs} loop
         enddo alpha1loop ! over alpha1={1,norbs} loop

         cumat_t(sigma1, sigma2, sigma3, sigma4) = ctmp
     enddo sigma4loop ! over sigma4={1,norbs} loop
     enddo sigma3loop ! over sigma3={1,norbs} loop
     enddo sigma2loop ! over sigma2={1,norbs} loop
     enddo sigma1loop ! over sigma1={1,norbs} loop

# if defined (Duliang)
     open(mytst, file='test-tumat.dat', status='unknown')
     do sigma1=1,norbs
     do sigma2=1,norbs
         do sigma3=1,norbs
         do sigma4=1,norbs
             if (abs(cumat_t(sigma1,sigma2,sigma3,sigma4)) .lt. eps6) cycle
             write(mytst, '(4i8, 2f17.10)') sigma1,sigma2,sigma3,sigma4, &
             cumat_t(sigma1, sigma2, sigma3, sigma4)
         enddo ! over gamma={1,norbs} loop
         enddo ! over delta={1,norbs} loop
     enddo ! over betta={1,norbs} loop
     enddo ! over alpha={1,norbs} loop
     close(mytst)
# endif /* Duliang */
     return
  end subroutine atomic_tran_cumat

  subroutine atomic_make_amtrx(norbs, cemat, somat, elocs, amtrx)
     use constants

     implicit none

! external arguments
! number of orbits
     integer, intent(in) :: norbs

! transformation matrix from {px, sz} to {j2, jz} basis
     complex(dp), intent(in) :: cemat(norbs, norbs)

! transformation matrix from {lz, sz} to {j2, jz} basis
     complex(dp), intent(in) :: somat(norbs, norbs)

! eigenvalue after diagonalizing the (cemat + somat)
     real(dp), intent(out) :: elocs(norbs)

! transformation matrix from {px, sz} to {j2, jz} basis
     complex(dp), intent(out) :: amtrx(norbs, norbs)

! auxiliary complex(dp) arrays
     complex(dp), allocatable :: ztmpa(:, :)

! local variables
     integer :: i, j

! allocate memory for auxiliary complex(dp) arrays
     allocate(ztmpa(norbs, norbs)); ztmpa = dcmplx(0.0D0, 0.0D0)

     amtrx = czero; ztmpa = cemat + somat
     if(sum(abs(ztmpa)) .lt. eps9)then
         do i = 1, norbs
             amtrx(i, i) = dcmplx(1.d0, 0.d0)
         enddo
     else
         call zmat_zheev(norbs, norbs, ztmpa, elocs, amtrx)
         ! added by sypeng, not general and should be modified!         
         amtrx = czero
         elocs = 0.d0
         do i = 1, norbs
             amtrx(i, i) = dcmplx(1.d0, 0.d0)
             elocs(i)    = dble(cemat(i, i))
         enddo
         ! sypeng end         
     endif

! deallocate memory for allocated auxiliary arrays
     if (allocated(ztmpa)) deallocate(ztmpa)

     return
  end subroutine atomic_make_amtrx
