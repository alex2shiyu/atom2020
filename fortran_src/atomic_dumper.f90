!=========================================================================!
! project : rambutan
! program : atomic_dumper
! history : May 07, 2011
! authors : duliang (duleung@gmail.com)
! purpose : this file is just for test purpose
! comment : 
!=========================================================================!
!>>> dump eigen-values of atomic Hamiltonian to file
  subroutine atomic_dump_heigs(ncfgs, heigs)
     use constants

     implicit none

! external arguments
! number of configurations
     integer :: ncfgs

! eigenvalues of atomic Hamiltonian
     real(dp) :: heigs(ncfgs)

! local variables
! loop index over basis
     integer :: ibas

! open output file
     open(mytmp, file="atom.eval.in", status="unknown")

! dumper data into output file
     do ibas=1,ncfgs
         write(mytmp, "(I8, F17.10)") ibas, heigs(ibas)
     enddo ! over ivec={1,ncfgs} loop

! close output file
     close(mytmp)

     return
  end subroutine atomic_dump_heigs

!>>> dump eigen-states of atomic Hamiltonian to file
  subroutine atomic_dump_heigv(norbs, ncfgs, heigv, invcd)
     use constants

     implicit none

! external arguments
! number of orbits
     integer, intent(in) :: norbs

! number of configurations
     integer, intent(in) :: ncfgs

! binary code 
     integer, intent(in) :: invcd(norbs, ncfgs)

! eigenvalues of "ntot" subspace
     complex(dp), intent(in) :: heigv(ncfgs, ncfgs)

! local variables
! loop index over basis
     integer :: ibas
     integer :: jbas

! open output file
     open(mytmp, file="atom.eigs.in", status="unknown")

! dumper data into output file
     do jbas=1,ncfgs
         do ibas=1,ncfgs
             if (abs(heigv(ibas, jbas)) .lt. eps6) cycle
             write(mytmp, "(2I8, 2F17.10, 4X, 14I1)") ibas, jbas, heigv(ibas, jbas), invcd(:, ibas)
         enddo ! over ibas={1,ncfgs} loop
     enddo ! over jbas={1,ncfgs} loop

! close output file
     close(mytmp)

     return
  end subroutine atomic_dump_heigv

  subroutine atomic_dump_jjmat(ncfgs, jjmat)
     use constants

     implicit none

! number of configurations
     integer, intent(in) :: ncfgs

! j2 operator matrix in fixed many body basis
     complex(dp), intent(in) :: jjmat(ncfgs, ncfgs)

! local variables
! loop index over basis
     integer :: ibas
     integer :: jbas

! open output file
     open(mytmp, file="solver.jjmat.out", status="unknown")

! dumper data into output file
     do jbas=1,ncfgs
         do ibas=1,ncfgs
             if (abs(jjmat(ibas, jbas)) .lt. eps6) cycle
             write(mytmp, "(2I8, 2F12.6)") ibas, jbas, jjmat(ibas, jbas)
         enddo ! over ibas={1,ncfgs} loop
     enddo ! over jbas={1,ncfgs} loop

! close output file
     close(mytmp)

      return
  end subroutine atomic_dump_jjmat

  subroutine atomic_dump_jzmat(ncfgs, jzmat)
     use constants

     implicit none

! number of configurations
     integer, intent(in) :: ncfgs

! jz operator matrix in fixed many body basis
     complex(dp), intent(in) :: jzmat(ncfgs, ncfgs)

! local variables
! loop index over basis
     integer :: ibas
     integer :: jbas

! open output file
     open(mytmp, file="solver.jzmat.out", status="unknown")

! dumper data into output file
     do jbas=1,ncfgs
         do ibas=1,ncfgs
             if (abs(jzmat(ibas, jbas)) .lt. eps6) cycle
             write(mytmp, "(2I8, 2F12.6)") ibas, jbas, jzmat(ibas, jbas)
         enddo ! over ibas={1,ncfgs} loop
     enddo ! over jbas={1,ncfgs} loop

! close output file
     close(mytmp)

      return
  end subroutine atomic_dump_jzmat

  subroutine atomic_dump_taeig(ncfgs, neigs, Heigs, jjeig, jzeig)
     use constants

     implicit none

! number of configurations
     integer, intent(in) :: ncfgs

! eigenvalue of jz operator matrix
     real(dp), intent(in) :: neigs(ncfgs)
     real(dp), intent(in) :: jjeig(ncfgs)
     real(dp), intent(in) :: jzeig(ncfgs)
     real(dp), intent(in) :: Heigs(ncfgs)

! local variables
! loop index over basis
     integer :: ibas

! open output file
     open(mytmp, file="solver.taeig.out", status="unknown")

! dumper data into output file
     write(mytmp, *) "ibas, heigs(ibas), neigs(ibas), jjeig(ibas), jzeig(ibas)"
     do ibas=1,ncfgs
         write(mytmp, "(I8, 4F17.10)") ibas, heigs(ibas), neigs(ibas), jjeig(ibas), jzeig(ibas)
     enddo ! over ibas={1,ncfgs} loop

! close output file
     close(mytmp)

      return
  end subroutine atomic_dump_taeig


!     call atomic_dump_taeig2(ncfgs, irrep_type, neigs, sgeigs, szeig, sseig, Heigs)
!  subroutine atomic_dump_taeig2(ncfgs, neigs, Heigs, jjeig, jzeig)
  subroutine atomic_dump_taeig2(ncfgs, irrep_type, neigs, sgeigs, szeig, sseig, Heigs)
     use constants

     implicit none

! number of configurations
     integer, intent(in) :: ncfgs

! eigenvalue of jz operator matrix
     integer , intent(in) :: irrep_type(ncfgs)
     real(dp), intent(in) :: neigs(ncfgs)
     real(dp), intent(in) :: sgeigs(ncfgs)
     real(dp), intent(in) :: szeig(ncfgs)
     real(dp), intent(in) :: sseig(ncfgs)
     real(dp), intent(in) :: Heigs(ncfgs)

! local variables
! loop index over basis
     integer :: ibas

! open output file
     open(mytmp, file="atom.goodquanum.dat", status="unknown")

! dumper data into output file
     write(mytmp, *) "ibas,    irrep_type,    neigs(ibas),    sgeigs(ibas),    sseig(ibas),    szeig(ibas),    Heigs(ibas)"
     do ibas=1,ncfgs
         write(mytmp, "(I8, I8, 5F17.10)") ibas, irrep_type(ibas), neigs(ibas), sgeigs(ibas), sseig(ibas), szeig(ibas), Heigs(ibas)
     enddo ! over ibas={1,ncfgs} loop

! close output file
     close(mytmp)

      return
  end subroutine atomic_dump_taeig2


! dump c4mat in many body space 
  subroutine atomic_dump_c4mat(ncfgs, c4mat)
     use constants

     implicit none

! number of configurations
     integer, intent(in) :: ncfgs

! l2 operator in many body basis
     complex(dp), intent(in) :: c4mat(ncfgs, ncfgs)

     open(mytmp, file="atom.c4mat.in")
         call zmat_dump2(mytmp, ncfgs, ncfgs, c4mat)
     close(mytmp)

     return
  end subroutine atomic_dump_c4mat
 
! dump c4mat in many body space 
  subroutine atomic_dump_c4mat2(ncfgs, c4mat, c4path)
     use constants

     implicit none

! number of configurations
     integer, intent(in) :: ncfgs

! l2 operator in many body basis
     complex(dp), intent(in) :: c4mat(ncfgs, ncfgs)

! file name
     character(len = *), intent(in) :: c4path

     open(mytmp, file= c4path)
         call zmat_dump2(mytmp, ncfgs, ncfgs, c4mat)
     close(mytmp)

     return
  end subroutine atomic_dump_c4mat2

  subroutine atomic_dump_irrep_flag(ncfgs, A, B, C,  c4path)
     use constants

     implicit none

! number of configurations
     integer, intent(in) :: ncfgs

! l2 operator in many body basis
     integer, intent(in) :: A(ncfgs)

! l2 operator in many body basis
     integer, intent(in) :: B(ncfgs)

! l2 operator in many body basis
     integer, intent(in) :: C(ncfgs)

! file name
     character(len = *), intent(in) :: c4path

! local variables
     integer    ::  i

     open(mytmp, file= c4path)
     do i = 1, ncfgs
         if(C(i) .eq. 1)then
             write(mytmp, '(I8,I8,A8)') A(i), B(i), 'A1'
         elseif(C(i) .eq. 2)then
             write(mytmp, '(I8,I8,A8)') A(i), B(i), 'A2'
         elseif(C(i) .eq. 3)then
             write(mytmp, '(I8,I8,A8)') A(i), B(i), 'B1'
         elseif(C(i) .eq. 4)then
             write(mytmp, '(I8,I8,A8)') A(i), B(i), 'B2'
         elseif(C(i) .eq. 5)then
             write(mytmp, '(I8,I8,A8)') A(i), B(i), 'E'
         else
             print *,"i, A(i), B(i), C(i)"
             print *,i, A(i), B(i), C(i)
             stop "something wrong happen in <atomic_dump_irrep_flag>!"
         endif
     enddo
     close(mytmp)

     return
  end subroutine atomic_dump_irrep_flag

  
  
! dump one dimension real vector in many body space 
  subroutine atomic_dump_rsgnmat(ncfgs, sgnmat, filepath)
     use constants

     implicit none

! number of configurations
     integer, intent(in) :: ncfgs

! l2 operator in many body basis
     real(8), intent(in) :: sgnmat(ncfgs)

! file name
     character(len = *), intent(in) :: filepath

! local variables
     integer :: i

     open(mytmp, file= filepath)
     do  i = 1, ncfgs
         write(mytmp, '(I8,F17.10)')i, sgnmat(i)
     enddo
     close(mytmp)

     return
  end subroutine atomic_dump_rsgnmat


! dump one dimension integer vector in many body space 
  subroutine atomic_dump_isgnmat(ncfgs, sgnmat, filepath)
     use constants

     implicit none

! number of configurations
     integer, intent(in) :: ncfgs

! l2 operator in many body basis
     integer, intent(in) :: sgnmat(ncfgs)

! file name
     character(len = *), intent(in) :: filepath

! local variables
     integer :: i

     open(mytmp, file= filepath)
     do  i = 1, ncfgs
         write(mytmp, '(I8,I8)')i, sgnmat(i)
     enddo
     close(mytmp)

     return
  end subroutine atomic_dump_isgnmat


! dump degeneracy information, not general.
  subroutine atomic_dump_deglist(ncfgs, num_tot, deg_list, filepath)
     use constants

     implicit none

! number of configurations
     integer, intent(in) :: ncfgs

! number of configurations
     integer, intent(in) :: num_tot

! l2 operator in many body basis
     integer, intent(in) :: deg_list(ncfgs)

! file name
     character(len = *), intent(in) :: filepath

! local variables
     integer :: i

     open(mytmp, file= filepath)
     write(mytmp, *)"the number of total vpms:   ", num_tot
     write(mytmp, *)"the number of every dimensionility of degeneracy subspace: " 
     do  i = 1, ncfgs
         if(deg_list(i) .ne. 0)then
             write(mytmp, '(I8,5x,I8)')i,deg_list(i) 
         endif
     enddo
     close(mytmp)

     return
  end subroutine atomic_dump_deglist



  subroutine atomic_dump_tamat(ncfgs, l2mat, s2mat, j2mat, lzmat, szmat, jzmat)
     use constants

     implicit none

! number of configurations
     integer, intent(in) :: ncfgs

! l2 operator in many body basis
     complex(dp), intent(in) :: l2mat(ncfgs, ncfgs)

! s2 operator in many body basis
     complex(dp), intent(in) :: s2mat(ncfgs, ncfgs)

! j2 operator in many body basis
     complex(dp), intent(in) :: j2mat(ncfgs, ncfgs)

! lz operator in many body basis
     complex(dp), intent(in) :: lzmat(ncfgs, ncfgs)

! sz operator in many body basis
     complex(dp), intent(in) :: szmat(ncfgs, ncfgs)

! jz operator in many body basis
     complex(dp), intent(in) :: jzmat(ncfgs, ncfgs)

     open(mytmp, file="atom.l2mat.in")
         call zmat_dump2(mytmp, ncfgs, ncfgs, l2mat)
     close(mytmp)

     open(mytmp, file="atom.s2mat.in")
         call zmat_dump2(mytmp, ncfgs, ncfgs, s2mat)
     close(mytmp)

     open(mytmp, file="atom.j2mat.in")
         call zmat_dump2(mytmp, ncfgs, ncfgs, j2mat)
     close(mytmp)

     open(mytmp, file="atom.lzmat.in")
         call zmat_dump2(mytmp, ncfgs, ncfgs, lzmat)
     close(mytmp)

     open(mytmp, file="atom.szmat.in")
         call zmat_dump2(mytmp, ncfgs, ncfgs, szmat)
     close(mytmp)

     open(mytmp, file="atom.jzmat.in")
         call zmat_dump2(mytmp, ncfgs, ncfgs, jzmat)
     close(mytmp)

     return
  end subroutine atomic_dump_tamat

  subroutine atomic_dump_vpmmat(ncfgs, vpmmat)
     use constants

     implicit none

! number of configurations
     integer, intent(in) :: ncfgs

! gutzwiller variational parameters matrix
     integer, intent(in) :: vpmmat(ncfgs, ncfgs)

! local variables
! loop index over basis
     integer :: ibas
     integer :: jbas

! open output file
     open(mytmp, file="atom.vpmsy.in", status="unknown")

! dumper data into output file
     do ibas=1,ncfgs
         do jbas=1,ncfgs
             if (vpmmat(ibas, jbas) .eq. 0) cycle
             write(mytmp, "(3I8)") ibas, jbas, vpmmat(ibas, jbas)
         enddo ! over ibas={1,ncfgs} loop
     enddo ! over jbas={1,ncfgs} loop

! close output file
     close(mytmp)

      return
  end subroutine atomic_dump_vpmmat

  subroutine atomic_dump_hmtrx(ncfgs, hmtrx)
     use constants

     implicit none

! external arguments
! number of configurations
     integer, intent(in) :: ncfgs

! atomic Hamiltonian matrix
     complex(dp), intent(in) :: hmtrx(ncfgs, ncfgs)

! local variables
! loop index over configurations
     integer :: ibas
     integer :: jbas

! open output file
     open(mytmp, file="atom.hmat.in", status="unknown")

! dumper data into output file
     do ibas=1,ncfgs
         do jbas=1,ncfgs
             if (abs(hmtrx(ibas, jbas)) .lt. eps6) cycle
             write(mytmp, "(2I8, 2F17.10)") ibas, jbas, Hmtrx(ibas, jbas)
         enddo ! over jbas={1,ncfgs} loop
     enddo ! over ibas={1,ncfgs} loop

! close output file
     close(mytmp)

     return
  end subroutine atomic_dump_hmtrx

  subroutine atomic_dump_amtrx(norbs, amtrx)
     use constants

     implicit none

! external arguments
! number of configurations
     integer, intent(in) :: norbs

! transformation matrix to make local Fock terms absent
     complex(dp), intent(in) :: amtrx(norbs, norbs)

! local variables
! loop index over configurations
     integer :: iorb
     integer :: jorb

! open output file
     open(mytmp, file="atom.amat.in", status="unknown")

! dumper data into output file
     do jorb=1,norbs
         do iorb=1,norbs
             if (abs(amtrx(iorb, jorb)) .lt. eps6) cycle
             write(mytmp, "(2I8, 2F20.14)") iorb, jorb, amtrx(iorb, jorb)
         enddo ! over jbas={1,ncfgs} loop
     enddo ! over ibas={1,ncfgs} loop

! close output file
     close(mytmp)

     return
  end subroutine atomic_dump_amtrx

  subroutine atomic_dump_elocs(norbs, elocs)
     use constants

     implicit none

! external arguments
! number of configurations
     integer, intent(in) :: norbs

! local onsite energy levels (should be diagonal)
     real(dp), intent(in) :: elocs(norbs)

! local variables
! loop index over configurations
     integer :: iorb

! open output file
     open(mytmp, file="atom.eloc.in", status="unknown")

! dumper data into output file
     do iorb=1,norbs
         write(mytmp, "(I8, F17.10)") iorb, elocs(iorb)
     enddo ! over ibas={1,ncfgs} loop

! close output file
     close(mytmp)

     return
  end subroutine atomic_dump_elocs

  subroutine atomic_dump_basis(norbs, ncfgs, basis, invsn, invcd)
     use constants

     implicit none

! external arguments
! number of orbits
     integer, intent(in) :: norbs

! number of configurations
     integer, intent(in) :: ncfgs

! decimal representation of Fock basis
     integer, intent(in) :: basis(ncfgs)

! serial number of decimal represented Fock basis
     integer, intent(in) :: invsn(0:2**norbs-1)

! binary representation of Fock basis
     integer, intent(in) :: invcd(norbs, ncfgs)

! local variables
! loop index over configurations
     integer :: ibas

! open output file
     open(mytmp, file="atom.basis.in", status="unknown")

! dumper data into output file
     do ibas=1,ncfgs
         write(mytmp, "(3I8,3X,14I2)") ibas, basis(ibas), &
         invsn(basis(ibas)), invcd(1:norbs, ibas)
     enddo ! over ibas={1,nconf(ntot)} loop

! close output file
     close(mytmp)

     return
  end subroutine atomic_dump_basis
