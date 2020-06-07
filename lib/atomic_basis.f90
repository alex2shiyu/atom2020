!========================================================================!
! project : clematis
! program : atomic_build_basis
! history : Apr 29, 2011
! authors : duliang (emial:duleung@gmail.com)
! purpose : build basis in "ntots" subspace or global space
! comment : 
!=========================================================================!
  subroutine atomic_make_basis(norbs, totncfgs, ncfgs, ntots, nmin, nmax, nstat, iprint, basis, invcd, invsn)

     implicit none

! external variables
! number of orbits
     integer, intent(in) :: norbs

! number of total configurations(totncfgs should be set 2**norbs as input)
! I dont't use 2**norbs in code instead of totncfgs because f2py will recognize "**" 
! as pointer by gcc.
     integer, intent(in) :: totncfgs

     ! number of Fock basis
     integer, intent(in) :: ncfgs

! number of total electrons
     integer, intent(in) :: ntots

! min number of electron occupancy
     integer, intent(in) :: nmin

! max number of electron occupancy
     integer, intent(in) :: nmax

! number of configuration
     integer, intent(in) :: nstat(0:norbs)

     integer, intent(in) :: iprint

! decimal representation of Fock basis
     integer, intent(out) :: basis(ncfgs)

! serial number of decimal represented Fock basis
     integer, intent(out) :: invsn(0:totncfgs-1)

! binary representation of Fock basis
     integer, intent(out) :: invcd(norbs, ncfgs)

! local variables
! loop index over decimal number
     integer :: inum

! loop index over fock basis
     integer :: ibas

! loop index over orbits
     integer :: ipos

! basis counter in global space
     integer :: icnt
     integer :: jcnt

! number of total electrons
     integer :: itot

! number of total electrons denoted by a decical number
     integer :: nbit

! auxiliary integer variables
     integer :: casea, caseb, casec
     integer :: ntiny, nlarg


!f2py intent(in) norbs 
!f2py intent(in) totncfgs
!f2py intent(in) ncfgs 
!f2py intent(in) ntots
!f2py intent(in) nmin
!f2py intent(in) nmax
!f2py intent(in) nstat
!f2py intent(in) iprint
!f2py intent(out) basis
!f2py intent(out) invsn
!f2py intent(out) invcd
!f2py depend(norbs) nstat
!f2py depend(ncfgs) basis
!f2py depend(totncfgs) invsn
!f2py depend(norbs,ncfgs) invcd


! initialize some variables
     icnt = 0; basis = 0
     jcnt = 0; invsn = 0; invcd = 0

! casea => construct basis sheet in only itot subspace (default)
! caseb => construct basis sheet in global state space (special)
!!!     casea = nstat(ntots); caseb = sum(nstat(0:norbs))
     casea = sum(nstat(0:ntots)); caseb = sum(nstat(0:norbs)); casec = sum(nstat(nmin:nmax))
     ! case I : many body basis of a subspace
     if (ncfgs .eq. casea) then
!!!         ntiny = ntots; nlarg = ntots
         ntiny = 0; nlarg = ntots
     endif ! back if (ncfgs .eq. casea) block

     ! case II: many body basis of global space
     if (ncfgs .eq. caseb) then
         ntiny = 0    ; nlarg = norbs
     endif ! back if (ncfgs .eq. caseb) block

     ! case III: many body basis of selective space
     if (ncfgs .eq. casec) then
         ntiny = nmin    ; nlarg = nmax
     endif ! back if (ncfgs .eq. caseb) block

     ! case III: severe error case
     if ((ncfgs.ne.casea) .and. (ncfgs.ne.caseb) .and. (ncfgs.ne.casec)) then
         stop "number of configurations ncfgs is wrong, modify dft.atom.in first"
     endif ! back if ((ncfgs.ne.casea) .and. (ncfgs.eq.caseb)) block

! construct basis(decimal) and "serial number(sn)" in Fock space
     itotloop: do itot=ntiny,nlarg

         do inum=0,2**norbs-1
             call verif( inum, norbs, nbit )
             if ( nbit .eq. itot ) then
                 icnt = icnt + 1; jcnt = jcnt + 1
                 basis(icnt) = inum; invsn(inum) = icnt
             endif ! back if ( nbit .eq. itot ) block
         enddo ! over inum={0,2**norbs-1} loop

!# very important, check number of configurations in itot subspace
         if (jcnt .ne. nstat(itot)) then
             stop "error happened in subroutine atomic_build_basis"
         endif ! back if ( icnt .ne. nstat(itot) ) block
         jcnt = 0 ! reinitialize the counter for next subspace
        
     enddo itotloop ! over itot={ntiny,nlarg} loop

! construct inverse binary code from a decimal number 
     do ibas=1,ncfgs
         do ipos=1,norbs
             if( btest(basis(ibas), ipos-1) ) invcd(ipos, ibas) = 1
         enddo ! over ipos={1,norbs} loop
     enddo ! over ibas={1,ncfgs} loop

! dump atomic configurations to file "atom.basis.in"
     if (iprint .ge. 1)then
         call atomic_dump_basis(norbs,totncfgs, ncfgs, basis, invsn, invcd)
     endif

     return
  end subroutine atomic_make_basis



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
!>>> electron number denoted by a decimal number
  subroutine verif( inum, norbs, nbits )
     implicit none

! external arguments
! decimal number
     integer, intent(in)  :: inum

! number of orbits
     integer, intent(in)  :: norbs

! number of total electrons
     integer, intent(out) :: nbits

! local variables
! local index over orbits
     integer :: ipos

! number of total electrons
     nbits = 0
     do ipos=1,norbs
         if( btest(inum, ipos-1) ) nbits = nbits + 1
     enddo ! over ipos={1,norbs} loop

     return
  end subroutine verif
