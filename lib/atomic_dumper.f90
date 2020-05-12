  subroutine atomic_dump_basis(norbs, totncfgs, ncfgs, basis, invsn, invcd)

     implicit none

! external arguments
! number of orbits
     integer, intent(in) :: norbs

! number of total configurations(totncfgs should be set 2**norbs as input)
! I dont't use 2**norbs in code instead of totncfgs because f2py will recognize "**" 
! as pointer by gcc.
     integer, intent(in) :: totncfgs

! number of configurations
     integer, intent(in) :: ncfgs

! decimal representation of Fock basis
     integer, intent(in) :: basis(ncfgs)

! serial number of decimal represented Fock basis
     integer, intent(in) :: invsn(0:totncfgs-1)

! binary representation of Fock basis
     integer, intent(in) :: invcd(norbs, ncfgs)

!f2py intent(in)  norbs
!f2py intent(in)  totncfgs
!f2py intent(in)  ncfgs
!f2py intent(in)  basis
!f2py intent(in)  invsn
!f2py intent(in)  invcd
!f2py depend(ncfgs)  basis
!f2py depend(norbs,ncfgs)  invcd
!f2py depend(totncfgs) invsn



! local variables
! loop index over configurations
     integer :: ibas

! open output file
     open(111, file="atom.basis.in", status="unknown")

! dumper data into output file
     do ibas=1,ncfgs
         write(111, "(3I8,3X,14I2)") ibas, basis(ibas), &
         invsn(basis(ibas)), invcd(1:norbs, ibas)
     enddo ! over ibas={1,nconf(ntot)} loop

! close output file
     close(111)

     return
  end subroutine atomic_dump_basis
