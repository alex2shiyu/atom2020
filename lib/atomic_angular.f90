! Operator transformation from single particle basis to many body basis
  subroutine atomic_make_sp2np(norbs, ncfgs, totncfgs, basis, invsn, invcd, zmat_s, zmat_n)
     implicit none

! number of orbits
     integer, intent(in) :: norbs

! number of configurations
     integer, intent(in) :: ncfgs

! number of total configurations of the whole shell
     integer, intent(in) :: totncfgs

! decimal represented basis
     integer, intent(in) :: basis(ncfgs)

! serial number of a decimal represented basis
     integer, intent(in) :: invsn(0:totncfgs-1)

! binary represented basis
     integer, intent(in) :: invcd(norbs, ncfgs)

! single particle operator in single particle basis
     complex(8), intent(in) :: zmat_s(norbs, norbs)

! single particle operator in many particle basis
     complex(8), intent(out) :: zmat_n(ncfgs, ncfgs)

! local variables
! loop index over orbits
     integer :: iorb
     integer :: alpha
     integer :: betta

! loop index over configurations
     integer :: ibas
     integer :: jbas

! sign change due to anti-commute relation
     integer :: isgn

! auxiliary integer variables
     integer :: jold, jnew
     integer :: code(norbs)

! auxiliary real(dp) variables
     real(8) :: dsgn

!f2py   intent(in)     norbs
!f2py   intent(in)     ncfgs
!f2py   intent(in)     totncfgs
!f2py   intent(in)     basis
!f2py   intent(in)     invsn
!f2py   intent(in)     invcd
!f2py   intent(in)     zmat_s
!f2py   intent(out)    zmat_n
!f2py   depend(ncfgs)  basis
!f2py   depend(totcfgs)  invsn
!f2py   depend(norbs,ncfgs)  invcd
!f2py   depend(norbs)  zmat_s
!f2py   depend(ncfgs)  zmat_n


! initialize zmat_n to be zero
     zmat_n = dcmplx(0.0d0, 0.0d0)

! main loop over many body basis
     jbasloop: do jbas=1,ncfgs
         alphaloop: do alpha=1,norbs
         bettaloop: do betta=1,norbs

             ! initialize some variables
             if (abs(zmat_s(alpha, betta)) .lt. 1.0D-10) cycle
             isgn = 0; jold = basis(jbas); code = invcd(:, jbas)

             ! simulate one eliminate operator
             if (code(betta) .eq. 1) then
                 do iorb=1,betta-1
                     isgn = isgn + code(iorb)
                 enddo ! over iorb={1,norbs} loop
                 code(betta) = 0

             ! simulate one construct operator
             if (code(alpha) .eq. 0) then
                 do iorb=1,alpha-1
                     isgn = isgn + code(iorb)
                 enddo ! over iorb={1,norbs} loop
                 code(alpha) = 1

                 ! determine the row number and matrix element
                 dsgn = 1.0d0; if (mod(isgn, 2) .eq. 1) dsgn = - 1.0d0
                 jnew = jold - 2**(betta-1) + 2**(alpha-1); ibas = invsn(jnew)

                 if (ibas .eq. 0) stop "severe error happened in atomic_make_sp2np"
                 zmat_n(ibas, jbas) = zmat_n(ibas, jbas) + zmat_s(alpha, betta) * dsgn

             endif ! back if (code(betta) .eq. 0) block
             endif ! back if (code(alpha) .eq. 1) block

         enddo bettaloop ! over betta={1,norbs} loop
         enddo alphaloop ! over alpha={1,norbs} loop
     enddo jbasloop ! over jbas={1,ncfgs} loop

     return
  end subroutine atomic_make_sp2np
