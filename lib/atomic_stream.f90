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
