!>>> create one electron on ipos of |jold) to deduce |jnew)
  subroutine atomic_construct(ipos, jold, jnew, isgn)

     implicit none

! external argument
! position number (serial number of orbit)
     integer, intent(in) :: ipos

! old Fock state and new Fock state
     integer, intent(in ):: jold
     integer, intent(out):: jnew

! sgn due to anti-commute relation between fernions
     integer, intent(out):: isgn

! local variables
! loop index over orbit
     integer :: iorb

     if (btest(jold, ipos-1) .eqv. .true.) then
         stop "severe error happened in atomic_construct"
     endif ! back if (btest(jold, ipos-1) .eqv. .true.) block

     isgn = 0
     do iorb=1,ipos-1
        if (btest(jold, iorb-1)) isgn = isgn + 1
     enddo ! over i={1,ipos-1} loop
     isgn = mod(isgn, 2)

     isgn = (-1)**isgn
     jnew = jold + 2**(ipos-1)

     return
  end subroutine atomic_construct

!>>> destroy one electron on ipos of |jold) to deduce |jnew)
  subroutine atomic_eliminate(ipos, jold, jnew, isgn)
     implicit none

! external argument
! position number (serial number of orbit)
     integer, intent(in)  :: ipos

! old Fock state and new Fock state
     integer, intent(in ) :: jold
     integer, intent(out) :: jnew

! sgn due to anti-commute relation between fernions
     integer, intent(out) :: isgn

! local variables
! loop index over orbit
     integer :: iorb

     if (btest(jold, ipos-1) .eqv. .false.) then
         stop "severe error happened in atomic_eliminate"
     endif ! back if (btest(jold, ipos-1) .eqv. .false.) block

     isgn = 0
     do iorb=1,ipos-1
         if (btest(jold, iorb-1)) isgn = isgn + 1
     enddo ! over i={1,ipos-1} loop
     isgn = mod(isgn, 2)

     isgn = (-1)**isgn
     jnew = jold - 2**(ipos-1)

     return
  end subroutine atomic_eliminate

!>>> build matrix for construct operator
  subroutine atomic_build_fpadd(norbs, ncfgs, basis, invsn, fpadd)
     use constants

     implicit none

! external arguments
! number of orbits
     integer, intent(in) :: norbs

! number of configurations
     integer, intent(in) :: ncfgs

! decimal represented Fock basis
     integer, intent(in) :: basis(ncfgs)

! inverse serial number of Fock basis
     integer, intent(in) :: invsn(0:2**norbs-1)

! matrix form of construct operator
     complex(dp), intent(out) :: fpadd(ncfgs, ncfgs, norbs)

! local variables
! loop index ober orbits
     integer :: iorb

! sign change due to commute relation
     integer :: isgn

! auxiliary integer variables
     integer :: jold, jnew

! loop index over configurations
     integer :: ibas, jbas

! initialize fpadd to be czero
     fpadd = czero

! this program works only for global space
     if ( 2**norbs .ne. ncfgs) then
         stop "subroutine atomic_build_fpadd works only for global space"
     endif ! back if (2**norbs .ne. ncfgs) block

! setup matrix for construct operator
     do iorb=1,norbs
         do jbas=1,ncfgs
             jold = basis(jbas)
             if (btest(jold, iorb-1) .eqv. .false.) then
                 call atomic_construct(iorb, jold, jnew, isgn)
                 ibas = invsn(jnew)
                 fpadd(ibas, jbas, iorb) = dble(isgn)
             endif
         enddo ! over jbas={1,ncfgs} loop
     enddo ! over iorb={1,norbs} loop

# if defined (duliang)
     open(mytst, file="solver.fpadd.out")
     call imat_dump3(mytst, ncfgs, ncfgs, norbs, fpadd)
     close(mytst)
# endif /* duliang */
     return
  end subroutine atomic_build_fpadd

!>>> build matrix form of eliminate operator
  subroutine atomic_build_fpcut(norbs, ncfgs, basis, invsn, fpcut)
     use constants

     implicit none

! external arguments
! number of orbits
     integer, intent(in) :: norbs

! number of configurations
     integer, intent(in) :: ncfgs

! decimal represented Fock basis
     integer, intent(in) :: basis(ncfgs)

! inverse serial number of Fock basis
     integer, intent(in) :: invsn(0:2**norbs-1)

! matrix form of eliminate operator
     complex(dp), intent(out) :: fpcut(ncfgs, ncfgs, norbs)

! local variables
! loop index ober orbits
     integer :: iorb

! sign change due to commute relation
     integer :: isgn

! auxiliary integer variables
     integer :: jold, jnew

! loop index over configurations
     integer :: ibas, jbas

! initialize fpcut to be czero
     fpcut = czero

! this program works only for global space
     if (2**norbs .ne. ncfgs) then
         stop "subroutine atomic_build_fpcut works only for global space"
     endif ! back if (2**norbs .ne. ncfgs) block

! setup matrix for eliminate operator
     do iorb=1,norbs
         do jbas=1,ncfgs
             jold = basis(jbas)
             if (btest(jold, iorb-1) .eqv. .true.) then
                 call atomic_eliminate(iorb, jold, jnew, isgn)
                 ibas = invsn(jnew)
                 fpcut(ibas, jbas, iorb) = dble(isgn)
             endif
         enddo ! over jbas={1,ncfgs} loop
     enddo ! over iorb={1,norbs} loop

# if defined (Duliang)
     open(mytst, file="solver.fpcut.out")
     call zmat_dump3(mytst, ncfgs, ncfgs, norbs, fpcut)
     close(mytst)
# endif /* Duliang */
     return
  end subroutine atomic_build_fpcut

  subroutine atomic_dump_fmtrx(norbs, ncfgs, eval, nval, fmat)
     use constants

     implicit none

! number of orbits
     integer, intent(in) :: norbs

! number of configurations
     integer, intent(in) :: ncfgs

     real(dp), intent(in) :: eval(ncfgs)
     real(dp), intent(in) :: nval(ncfgs)

! matrix form of construct operator
     complex(dp), intent(in) :: fmat(ncfgs, ncfgs, norbs)

! local variables
     integer :: i, j, k

     real(dp) :: sxyz 
     sxyz = 0.0d0

! open data file: atom.cix
     open(mytmp, file='atom.cix', form='formatted', status='unknown')

! write eigenvalues
     write(mytmp,'(a)') '# eigenvalues: index | energy | occupy | spin'
     do i=1,ncfgs
         write(mytmp,'(i5,3f16.8)') i, eval(i), nval(i), sxyz
     enddo ! over i={1,ncfgs} loop

! write f matrix element
     write(mytmp,'(a)') '# f matrix element: alpha | beta | orbital | fmat'
     do i=1,norbs
         do j=1,ncfgs
             do k=1,ncfgs
                 write(mytmp,'(3i5,f16.8)') k, j, i, real(fmat(k,j,i))
             enddo ! over k={1,ncfgs} loop
         enddo ! over j={1,ncfgs} loop
     enddo ! over i={1,norbs} loop

! close input file
     close(mytmp)

     return
  end subroutine atomic_dump_fmtrx
