!=========================================================================!
! project : strawberry
! program : atomic_build_Hmtrx
! history : 09/29/2011
! author  : xidai and duliang (email:duleung@gmail.com)
! purpose : construct atomic hamiltonian matrix
! comment : 
!=========================================================================!
  subroutine atomic_make_hmtrx(norbs, totncfgs, ncfgs, state, invcd, invsn, eimp, umat, iprint, hmat)

     implicit none

! external arguments
! number of orbits
     integer, intent(in) :: norbs

! number of Fock basis
     integer, intent(in) :: totncfgs

! number of Fock basis
     integer, intent(in) :: ncfgs

! decimal representation of Fock basis
     integer, intent(in) :: state(ncfgs)

! serial number of decimal represented Fokc basis
     integer, intent(in) :: invsn(0:totncfgs-1)

! binary representation of Fock basis
     integer, intent(in) :: invcd(norbs, ncfgs)

! impurity energy level
     complex(kind=8), intent(in) :: eimp(norbs, norbs)

! general coulomb interaction matrix
     complex(kind=8), intent(in) :: umat(norbs, norbs, norbs, norbs)

     integer, intent(in) :: iprint

! Hamiltonian matrix
     complex(kind=8), intent(out) :: hmat(ncfgs, ncfgs)

! local variables
! loop index
     integer :: i

! loop index over orbits
     integer :: iorb
     integer :: jorb

! sign change due to fermion anti-commute relation
     integer :: isgn

! loop index over Fock basis
     integer :: ibas, jbas

! new basis state after four fermion operation
     integer :: knew

! index for general interaction matrix
     integer :: alpha, betta
     integer :: delta, gamma

! binary code representation of a state
     integer :: code(norbs)

! auxiliary complex(dp) variables
     complex(kind=8) :: ztmpa

!f2py intent(in)  norbs
!f2py intent(in)  totncfgs
!f2py intent(in)  ncfgs
!f2py intent(in)  state
!f2py intent(in)  invsn
!f2py intent(in)  invcd
!f2py intent(in)  eimp
!f2py intent(in)  umat
!f2py intent(in)  iprint
!f2py intent(out) hmat
!f2py depend(ncfgs)       state
!f2py depend(totncfgs)    invsn
!f2py depend(norbs,ncfgs) invcd
!f2py depend(norbs)       eimp
!f2py depend(norbs)       umat
!f2py depend(ncfgs)       hmat



! initialize hmat
     hmat(1:ncfgs, 1:ncfgs) = dcmplx(0.d0,0.d0)

! local energy level Hamitonian
!     do ibas=1,ncfgs
!        code(1:norbs) = invcd(1:norbs, ibas)
!
!        do iorb=1,norbs
!           hmat(ibas, ibas) = hmat(ibas, ibas) + eloc(iorb) * code(iorb)
!        enddo ! over iorb={1,norbs} loop
!
!     enddo ! over ibas={1,ncfgs} loop

!=========================================================================!
! two fermion operator terms (crystalline electric field)
!=========================================================================!
     do jbas=1,ncfgs
         alploop: do alpha=1,norbs
         betloop: do betta=1,norbs

             isgn = 0
             knew = state(jbas)
             code(1:norbs) = invcd(1:norbs, jbas)
             if ( abs(eimp(alpha, betta)) .lt. 1E-8 ) cycle

! simulate one eliminate operator
             if (code(betta) == 1) then
                 do i=1,betta-1
                     if (code(i) == 1) isgn = isgn + 1
                 enddo 
                 code(betta) = 0

! simulate one construct operator
                 if (code(alpha) == 0) then
                     do i=1,alpha-1
                         if (code(i) == 1) isgn = isgn + 1
                     enddo
                     code(alpha) = 1

! determine the column number and hamiltonian matrix elememt
                     knew = knew - 2**(betta-1)
                     knew = knew + 2**(alpha-1)
 
                     isgn  = mod(isgn, 2)
                     ibas = invsn(knew)
                     if (ibas == 0) stop "error while determining row1"
                     hmat(ibas,jbas) = hmat(ibas,jbas) + eimp(alpha,betta) * (-1.0d0)**isgn 

                 endif ! back if (code(alpha) == 0) block
             endif ! back if (betta == 1) block

         enddo betloop ! over betta={1,norbs} loop
         enddo alploop ! over alpha={1,norbs} loop
     enddo ! over jbas={1,nbas} loop


! four fermion terms in local Hamiltonian
     do jbas=1,ncfgs

        alphaloop : do alpha=1,norbs
        bettaloop : do betta=1,norbs
        gammaloop : do gamma=1,norbs
        deltaloop : do delta=1,norbs

            isgn = 0
            knew = state(jbas)
            code(1:norbs) = invcd(1:norbs, jbas)
            !# very important if single particle basis rotated
            if ((alpha .eq. betta) .or. (delta .eq. gamma)) cycle
            if ( abs(umat(alpha,betta,delta,gamma)) .lt. 1E-8 ) cycle

! simulate two eliminate operator
            if ((code(delta) == 1) .and. (code(gamma) == 1)) then
                do i=1,gamma-1
                    if(code(i) == 1) isgn = isgn + 1
                enddo ! over i={1,gamma-1} loop
                code(gamma) = 0

                do i=1,delta-1
                    if(code(i) == 1) isgn = isgn + 1
                enddo ! over i={1,delta-1} loop
                code(delta) = 0

! simulate two construct operator
                if ((code(alpha) == 0) .and. (code(betta) == 0)) then
                    do i=1,betta-1
                        if(code(i) == 1) isgn = isgn + 1
                    enddo ! over i={1,betta-1} loop
                    code(betta) = 1

                    do i=1,alpha-1
                        if(code(i) == 1) isgn = isgn + 1
                    enddo ! over i={1,alpha-1} loop
                    code(alpha) = 1

! determine the column number and hamiltonian matrix elememt
                    knew = knew - 2**(gamma-1) - 2**(delta-1)
                    knew = knew + 2**(betta-1) + 2**(alpha-1)

                    ibas = invsn(knew)
                    isgn = mod(isgn, 2)
                    hmat(ibas,jbas) = hmat(ibas,jbas) + umat(alpha,betta,delta,gamma) * (-1.0d0)**isgn

! simplly check the fermion anti-commute relation
                    if ( isgn /= 0 ) then
                    !$  stop "something wrong in atomic_make_hmat, pls check carefull"
                    endif ! back if ( sgn /= 0 ) block

                endif ! back if ((code(delta) == 1) .and. (code(gamma) == 1)) block
            endif ! back if ((code(alpha) == 0) .and. (code(betta) == 0)) block

        enddo deltaloop ! over delta={gamma+1,norbs} loop
        enddo gammaloop ! over gamma={1,norbs-1} loop
        enddo bettaloop ! over betta={alpha+1,norbs} loop
        enddo alphaloop ! over alpha={1,norbs-1} loop

     enddo ! over jbas={1,ncfgs} loop

     if (iprint .eq. 3)then
         open(112, file='test-hamilt.dat', status='unknown')
         do jbas=1,ncfgs
             do ibas=1,ncfgs
                 if (abs(hmat(ibas, jbas)) .lt. 1E-6) cycle
                 write(112, '(2i8, 2f17.10)') ibas, jbas, hmat(ibas, jbas)
             enddo ! over ibas={1,ncfgs} loop
         enddo ! over jbas={1,ncfgs} loop
         close(112)
!        call atomic_dump_hmtrx(ncfgs, hmat)
     endif
     return
  end subroutine atomic_make_hmtrx
