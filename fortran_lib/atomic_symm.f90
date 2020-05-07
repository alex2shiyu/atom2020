  subroutine atomic_make_vpmmat(ncfgs, nneig, j2eig, jzeig, heigs, vpmmat)
     use constants

     implicit none

! external arguments
! number of configurations
     integer, intent(in) :: ncfgs

! eigenvalue of N operator
     real(dp), intent(in) :: nneig(ncfgs)

! eigenvalue of j2 operator
     real(dp), intent(in) :: j2eig(ncfgs)

! eigenvalue of jz operator
     real(dp), intent(in) :: jzeig(ncfgs)

! eigenvalue of H operator
     real(dp), intent(in) :: heigs(ncfgs)

! variational parameters matrix
     integer, intent(out) :: vpmmat(ncfgs, ncfgs)

! local variables
! loop index over configurations
     integer :: ibas
     integer :: jbas

! loop index over degenerancy subspace
     integer :: ideg

! number of degeneracy with 2*J+1
     integer :: ndegs
     integer :: mdegs

! counter of gutzwiller variational parameters
     integer :: ivpm

! auxiliary real(dp) variables
     real(dp) :: dtmpa
     real(dp) :: dtmpb

     real(8) :: jjeig(ncfgs)

! initialize vpmmat to be 0
     vpmmat = 0

     do ibas=1,ncfgs
         dtmpa = dble(j2eig(ibas))
         jjeig(ibas) = dsqrt(1.0D0+4.0D0*dtmpa)/2.0D0 - 0.5D0
     enddo ! over ibas={1,ncfgs}
!------------------------------------------------------------

!------------------------------------------------------------
! check the marshalling sequence of atomic eigenstates
!     ibas = 1
!     do while ( .true. )
!         ndegs = nint(dtwo * jjeig(ibas) + done)
!
!         do ideg=0,ndegs-1
!             dtmpa = jjeig(ibas+ideg) - jzeig(ibas+ideg)
!             if (abs(dtmpa - dble(ideg)) .gt. eps9) then
!                 write(60,*)'ibase,ideg,jjeig,jzeig'
!                 write(60,*)ibas,ideg,jjeig(ibas+ideg),jzeig(ibas+ideg)
!                 stop " severe error in atomic_make_vpmmat1"
!             endif
!         enddo ! over ideg={1,ndegs} loop
!
!         ibas = ibas + ndegs; if (ibas .gt. ncfgs) exit
!
!     enddo ! over while ( .true. ) loop

! setup variational parameters (diagonal parts)
     ivpm = 0; dtmpa = 1000.d0

!     do ibas=1,ncfgs
!         if (abs(heigs(ibas)-dtmpa) .gt. 1.0d-12) then
!             dtmpa = heigs(ibas); dtmpb = jjeig(ibas)
!             ivpm = ivpm + 1; vpmmat(ibas, ibas) = ivpm
!         else
!             if (abs(jjeig(ibas)-dtmpb) .lt. 1.0d-4) then
!                 vpmmat(ibas, ibas) = vpmmat(ibas-1, ibas-1)
!             else
!                 dtmpa = heigs(ibas); dtmpb = jjeig(ibas)
!                 ivpm = ivpm + 1; vpmmat(ibas, ibas) = ivpm
!             endif
!             if (ibas .eq. 1) stop "severe error in atomic_make_vpmmat2"
!         endif
!     enddo ! over ibas={1,ncfgs} loop
!
!! setup variational parameters (off-diagonal parts)
!     ibas = 1
!     do while(ibas .le. ncfgs)
!         ndegs = nint(dtwo*jjeig(ibas) + done)
!         jbas = ibas + ndegs
!
!         do while (jbas .le. ncfgs)
!             mdegs = nint(dtwo*jjeig(jbas) + done)
!
!             if (abs(nneig(jbas) - nneig(ibas)) .lt. eps9) then
!                 if (abs(jjeig(jbas) - jjeig(ibas)) .lt. eps9) then
!                     ivpm = ivpm + 1
!                     do ideg=1,ndegs
!                         vpmmat(ibas+ideg-1, jbas+ideg-1) = ivpm
!                         vpmmat(jbas+ideg-1, ibas+ideg-1) = ivpm
!                     enddo ! over ideg={1,ndegs} loop
!                 endif
!             endif
!
!             jbas = jbas + mdegs
!         enddo ! over jbas={1,ncfgs} loop
!         ibas = ibas + ndegs
!
!     enddo ! over ibas={1,ncfgs} loop

     write(*, "(A, I8)") "number of vpm is: ", ivpm
! check the Hermite symmetry of Project operator
     do ibas=1,ncfgs
         do jbas=1,ncfgs
             if (ibas .lt. jbas) cycle
             if (vpmmat(ibas, jbas) .ne. vpmmat(jbas, ibas)) then
                 stop "No Hermite symmetry of vpmmat in atomic_make_vpmmat"
             endif
         enddo ! over jbas={1,ncfgs} loop
     enddo ! over ibas={1,ncfgs} loop

! dump variational parameters into file
     call atomic_dump_vpmmat(ncfgs, vpmmat)

     return
  end subroutine atomic_make_vpmmat
  


  subroutine atomic_make_vpmmat2(ncfgs, Hdegs, nvpms, vpmmat)
     use constants

     implicit none

! external arguments
! number of configurations
     integer,  intent(in)  :: ncfgs

! degeneracy of atomic eigenvalues
     integer, intent(in)  :: hdegs(ncfgs)

! the number of vpms
     integer,  intent(out) :: nvpms

! variational parameters matrix
     integer,  intent(out) :: vpmmat(ncfgs, ncfgs)

! local variables
! loop index over configurations
     integer :: ibas
     integer :: jbas
     integer :: i, j

! loop index over degenerancy subspace
     integer :: ideg

! number of degeneracy 
     integer :: ndegs
     integer :: mdegs

! counter of gutzwiller variational parameters
     integer :: ivpm

! auxiliary real(dp) variables
     real(dp) :: dtmpa
     real(dp) :: dtmpb

     real(8)  :: jjeig(ncfgs)

     integer  :: ibase

! initialize vpmmat to be 0
     vpmmat = 0

! testing
     write(111,*)'hdegs',hdegs

! key part
     ibase = 0
     ibas  = 1
     ivpm  = 0
     do while ( ibas .le. ncfgs)
         ndegs = hdegs(ibas)
         do i = 1, ndegs
             do j = i, ndegs
                 if( j .eq. i)then
                     ivpm = ivpm + 1
                     vpmmat(j+ibas-1,i+ibas-1) = ivpm
                 else
                     ivpm = ivpm + 1
                     vpmmat(i+ibas-1,j+ibas-1) = ivpm
                     ivpm = ivpm + 1
                     vpmmat(j+ibas-1,i+ibas-1) = ivpm
                 endif
             enddo
         enddo
         ibas = ibas + hdegs(ibas)
     enddo

     nvpms = ivpm

! dump variational parameters into file
     call atomic_dump_vpmmat(ncfgs, vpmmat)

     return
  end subroutine atomic_make_vpmmat2
