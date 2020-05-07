!-------------------------------------------------------------------------!
! project : litchi
! program : atomic_clebsch_gordon
! history : oct 27, 2010
! authors : duliang (duleung@gmail.com)
! purpose : build Clebsch-Gordan coefficents matrix
! comment : ref: mazhongqi {wulixuezhongdequnlun} p223, 4.152
!-------------------------------------------------------------------------!
  real(8) function clebsch_gordan(j1, j2, j3, m1, m2, m3)
     implicit none

! orbital quantum number
     real(8), intent(in) :: j1, j2, j3

! magnetic quantum number
     real(8), intent(in) :: m1, m2, m3

! define some important parameters
! maxium integer number for factorial use
     integer, parameter :: maxnum = 40

! maxium loop for the summing process
     integer, parameter :: maxlop = 20

! array to record factorials => i!
     real(8) :: facts(0:maxnum)

! array to record phase (-1)^i
     real(8) :: phase(0:maxnum)

! loop index
     integer :: i, j

! temperary variables to record numerator and denominator of a fraction
     real(8) :: numer, denom

! prepare some utility arrays
     facts(0) = 1.0d0
     phase(0) = 1.0d0

     do i=1,maxnum
         facts(i) = facts(i-1) * dble( i)
         phase(i) = phase(i-1) * dble(-1)
     enddo ! over i={1,maxnum} loop

     clebsch_gordan = 0.0d0

     do i=0,maxlop
         if ( j2 + m1 - j3 + real(i) < 0.0d0 ) cycle
         if ( j1 - m1 - real(i) < 0.0d0 ) cycle
         if ( j3 - m3 - real(i) < 0.0d0 ) cycle
         numer = facts(nint(j3+j2-m1-real(i))) * facts(nint(j1+m1+real(i))) &
               * phase(nint(j1-m1+real(i)))
         denom = facts(nint(real(i))) * facts(nint(j1-m1-real(i)))          &
               * facts(nint(j3-m3-real(i))) * facts(nint(m1+j2-j3+real(i)))
         
         clebsch_gordan = clebsch_gordan + (numer / denom)
     enddo ! over i={0,maxlop} loop

     numer = facts(nint(j3+m3)) * facts(nint(j3-m3))         & 
           * facts(nint(j1-m1)) * facts(nint(j2-m2))         &
           * facts(nint(j1+j2-j3)) * (2.d0*j3+1.d0)

     denom = facts(nint(j1+m1)) * facts(nint(j2+m2))         &
            * facts(nint(j3-j1+j2)) * facts(nint(j3+j1-j2))  &
            * facts(nint(j1+j2+j3+1.d0))

     clebsch_gordan = clebsch_gordan * dsqrt(numer / denom)

     return
  end function clebsch_gordan
