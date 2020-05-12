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
subroutine atomic_make_cumat(norbs, Uc, Uv, Jz, Js, Jp, cumat )
     implicit none
! number of orbits
     integer, intent(in) :: norbs
! general coulomb interaction parameters
     real(kind=8), intent(in) :: Uc
     real(kind=8), intent(in) :: Uv
     real(kind=8), intent(in) :: Jz
     real(kind=8), intent(in) :: Js
     real(kind=8), intent(in) :: Jp
! general coulomb interaction matrix
     complex(kind=8), intent(out) :: cumat(norbs, norbs, norbs, norbs)

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
     real(kind=8) :: res

!f2py intent(in) norbs 
!f2py intent(in) Uc 
!f2py intent(in) Uv 
!f2py intent(in) Jz 
!f2py intent(in) Js 
!f2py intent(in) Jp
!f2py intent(out) cumat
!f2py depend(norbs) cumat

! initialize cumat to zero
     cumat = dcmplx(0.d0,0.d0)

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
            res = 0.d0

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

     open(11, file='test-cumat.dat', status='unknown')
     do alpha=1,norbs
     do betta=1,norbs
         do delta=1,norbs
         do gamma=1,norbs
             if ( abs(cumat(alpha, betta, delta, gamma)) .lt. 1e-6 ) cycle
             write(11, '(4i8, 2f17.10)') alpha, betta, delta, gamma, &
             cumat(alpha, betta, delta, gamma)
         enddo ! over gamma={1,norbs} loop
         enddo ! over delta={1,norbs} loop
     enddo ! over betta={1,norbs} loop
     enddo ! over alpha={1,norbs} loop
     close(11)
     return
end subroutine atomic_make_cumat


  subroutine atomic_tran_cumat(norbs, amtrx, cumat, cumat_t)

     implicit none

! number of orbits
     integer, intent(in) :: norbs

! transformation matrix from orginal basis to natural basis
     complex(kind=8), intent(in) :: amtrx(norbs, norbs)

! coefficents matrix for generalized interaction U in orginal basis
     complex(kind=8), intent(in) :: cumat(norbs, norbs, norbs, norbs)
 
! coefficents matrix for generalized interaction U in natural basis
     complex(kind=8), intent(out) :: cumat_t(norbs, norbs, norbs, norbs)

! local varoables
! loop index over orbits in orginal single particle basis
     integer :: alpha1, alpha2
     integer :: alpha3, alpha4

! loop index over orbits in natural single particle basis
     integer :: sigma1, sigma2
     integer :: sigma3, sigma4

! auxiliary complex(dp) variables
     complex(kind=8) :: ctmp

!f2py intent(in)  norbs
!f2py intent(in)  amtrx
!f2py intent(in)  cumat
!f2py intent(out) cumat_t
!f2py depend(norbs) amtrx
!f2py depend(norbs) cumat
!f2py depend(norbs) cumat_t


! initialize cumat_t to be zero
     cumat_t = dcmplx(0.0D0, 0.0D0)

     sigma1loop: do sigma1=1,norbs
     sigma2loop: do sigma2=1,norbs
     sigma3loop: do sigma3=1,norbs
     sigma4loop: do sigma4=1,norbs
         ctmp = dcmplx(0.0D0, 0.0D0)

         alpha1loop: do alpha1=1,norbs
         alpha2loop: do alpha2=1,norbs
         alpha3loop: do alpha3=1,norbs
         alpha4loop: do alpha4=1,norbs
             if (abs(cumat(alpha1, alpha2, alpha3, alpha4)) .lt. 1E-8) cycle
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

!!# if defined (duliang)
!!     open(mytmp, file='solver.tumat.out')
!!     call zmat_dump4(mytmp, norbs, norbs, norbs, norbs, cumat_t)
!!     close(mytmp)
!!# endif /* duliang */
     return
  end subroutine atomic_tran_cumat

!  subroutine atomic_state(norbs,nstat)
!
!     implicit none
!! number of orbitals
!     integer, intent(in)  :: norbs
!
!! number of configuration in each subspace
!     integer, intent(out) :: nstat(0:norbs)
!
!! local variables
!! loop index over good quantum number
!     integer :: ibit
!! loop index over orbits
!     integer :: iorb
!     integer :: jorb
!
!     do ibit=0,norbs
!         nstat(ibit) = state_pick(ibit, norbs)
!     enddo ! over ibit={0,norbs} loop 
!  end subroutine atomic_state
!
!  
!  !>>> calculate combination algebra 
!  function state_pick(ntiny, nlarg) result(value)
!     implicit none
!
!! external variables
!     integer, intent(in) :: ntiny
!     integer, intent(in) :: nlarg
!
!! local variables
!     integer :: i
!
!! auxiliary integer variable
!     integer :: nlow
!
!! numberator of the combination algebra
!     real(8) :: numer
!
!! denominator of the combination algebra
!     real(8) :: denom
!
!! result value of the combination algebra
!     integer :: value
!
!! transform the combination algebra
!     nlow = min(ntiny, nlarg-ntiny)
!
!! numerator in combination algebra
!     numer = 1.0D0
!     do i=nlarg-nlow+1,nlarg
!        numer = numer * dble(i)
!     enddo ! over i={nlarg-nlow+1,nlarg} loop
!
!! denominator in combination algebra
!     denom = 1.0D0
!     do i=1,nlow
!        denom = denom * dble(i)
!     enddo ! over i={1,nlow} loop
!
!! result value
!     value = nint(numer / denom)
!
!     return
!  end function state_pick

