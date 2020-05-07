  subroutine atomic_make_phase(ncfgs, jxmat, jymat, heigv)
     use constants

     implicit none

! external arguments
! number of configurations
     integer, intent(in) :: ncfgs

! eigenvalue of j2 operator
     complex(dp), intent(in) :: jxmat(ncfgs, ncfgs)

! eigenvalue of jz operator
     complex(dp), intent(in) :: jymat(ncfgs, ncfgs)

! eigenvector of atomic Hamiltonian matrix
     complex(dp), intent(inout) :: heigv(ncfgs, ncfgs)
 
! variational parameters matrix
     real(dp) :: pfact(ncfgs)

! local variables
! loop index over configurations
     integer :: ibas
     integer :: jbas

! loop index over degenerancy subspace
     integer :: ideg

! number of degeneracy with 2*J+1
     integer :: ndegs

! auxiliary real(dp) variables
     real(dp) :: dtmpa
     real(dp) :: dtmpb
     real(dp) :: dsgn(ncfgs)

! auxiliary complex(dp) matrix
     complex(dp) :: zmata(ncfgs, ncfgs)
     complex(dp) :: zmatb(ncfgs, ncfgs)
     complex(dp) :: jmmat(ncfgs, ncfgs)

! initialize 
     dtmpa = 0.d0
! initialize pfact to be positive 1
     pfact = 1.0d0

     jmmat = jxmat - dcmplx(0.0D0, 1.0D0) * jymat
     call zmat_zgemm0(ncfgs, jmmat, heigv, zmatb)
     call zmat_zgemm2(ncfgs, heigv, zmatb, zmata)
     call zmat_dump2(11, ncfgs, ncfgs, zmata)
!     do jbas=1,ncfgs
!         do ibas=1,ncfgs
!             if ((ibas-1) .eq. jbas) cycle
!             dtmpa = dtmpa + abs(zmata(ibas, jbas))
!         enddo ! over ibas={1,ncfgs} loop
!     enddo ! over jbas={1,ncfgs} loop
!     if (dtmpa .gt. 1.0D-8) stop "wrong structure of jmmat in atomic eignebasis"

     do jbas=1,ncfgs-1
         ibas = jbas + 1
         if ( abs(zmata(ibas, jbas)) .lt.  1.0D-8) pfact(ibas) =   1.0D0
         if (real(zmata(ibas, jbas)) .gt.  1.0D-8) pfact(ibas) =   pfact(jbas)
         if (real(zmata(ibas, jbas)) .lt. -1.0D-8) pfact(ibas) = - pfact(jbas)
     enddo ! over jbas={1,ncfgs} loop

     do jbas=1,ncfgs
         do ibas=1,ncfgs
             Heigv(ibas, jbas) = Heigv(ibas, jbas) * pfact(jbas)
         enddo ! over ibas={1,ncfgs} loop
     enddo ! over jbas={1,ncfgs} loop

# if defined (debug)
     call zmat_zgemm0(ncfgs, jmmat, heigv, zmatb)
     call zmat_zgemm2(ncfgs, heigv, zmatb, zmata)
     open(mytmp, file="pfact-debug.dat")
     do jbas=1,ncfgs
         do ibas=1,ncfgs
             if (abs(zmata(ibas, jbas)) .lt. 1.0D-8) cycle
             write(mytmp, "(2I8, 2F17.10)") ibas, jbas, zmata(ibas, jbas)
         enddo ! over ibas={1,ncfgs} loop
     enddo ! over jbas={1,ncfgs} loop
     close(mytmp)
# endif /* debug */
     return
  end subroutine atomic_make_phase
