! make the decimal number of Fock states 
subroutine MakeCg(norb,nmin,nmax,ncfgs,Cg)
    implicit none
    ! number of orb(with spin)
    integer, intent(in)    :: norb

    ! minimum occupancy
    integer, intent(in)    :: nmin

    ! maximum occupancy
    integer, intent(in)    :: nmax

    ! number of configurations
    integer, intent(in)    :: ncfgs

    ! decimal number
    integer, intent(out)   :: Cg(ncfgs)

    ! local variables
    integer   ::  icnt, jcnt, itot, inum, num_bit

!f2py   intent(in)     norb    
!f2py   intent(in)     nmin   
!f2py   intent(in)     nmax    
!f2py   intent(in)     ncfgs    
!f2py   intent(out)    Cg
!f2py   depend(ncfgs)  Cg


    icnt = 0
    jcnt = 0
    do itot=nmin,nmax
        do inum=0,2**norb-1
            call verif( inum, norb, num_bit )
            if ( num_bit .eq. itot ) then
                icnt = icnt + 1
                jcnt = jcnt + 1
                Cg(icnt) = inum
!                if(myid .eq. 0)write(19,*)'Cg',icnt,inum
            endif
        enddo
    enddo

end subroutine MakeCg


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

!f2py   intent(in)    inum
!f2py   intent(in)    norbs
!f2py   intent(out)   nbits

! number of total electrons
     nbits = 0
     do ipos=1,norbs
         if( btest(inum, ipos-1) ) nbits = nbits + 1
     enddo ! over ipos={1,norbs} loop

     return
  end subroutine verif
  
  ! a^b=c
  subroutine PowerInt(a,b,c)
      implicit none

      integer, intent(in)   ::  a

      integer, intent(in)   ::  b

      integer, intent(out)   ::  c

!f2py  intent(in)    a
!f2py  intent(in)    b
!f2py  intent(out)   c


      integer :: i
      
      c=1

      do i = 1, b
          c=c*a
      enddo

  end subroutine
