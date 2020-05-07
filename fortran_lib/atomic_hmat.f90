!=========================================================================!
! project : rambutan
! program : atomic_make_Hmtrx
! history : 12/05/2011
! author  : duliang (email:duleung@gmail.com)
! purpose : construct atomic hamiltonian matrix
! comment : 
!=========================================================================!
  subroutine atomic_make_Hmtrx(norbs, ncfgs, state, invcd, invsn, elocs, umat, hmat)
     use constants

     implicit none

! external arguments
! number of orbits
     integer, intent(in) :: norbs

! number of Fock basis
     integer, intent(in) :: ncfgs

! decimal representation of Fock basis
     integer, intent(in) :: state(ncfgs)

! serial number of decimal represented Fokc basis
     integer, intent(in) :: invsn(0:2**norbs-1)

! binary representation of Fock basis
     integer, intent(in) :: invcd(norbs, ncfgs)

! impurity energy level
     real(dp), intent(in) :: elocs(norbs)

! general coulomb interaction matrix
     complex(dp), intent(in) :: umat(norbs, norbs, norbs, norbs)

! Hamiltonian matrix
     complex(dp), intent(out) :: hmat(ncfgs, ncfgs)

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

! initialize hmat
     hmat(1:ncfgs, 1:ncfgs) = czero

! onsite energy terms
     do ibas=1,ncfgs
        do iorb=1,norbs
           if (invcd(iorb, ibas) .eq. 0) cycle
           Hmat(ibas, ibas) = Hmat(ibas, ibas) + elocs(iorb)
        enddo ! over iorb={1,norbs} loop
     enddo ! over ibas={1,ncfgs} loop

! four fermion term
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
            if ( abs(umat(alpha,betta,delta,gamma)) .lt. eps6 ) cycle

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
                    !!$  stop "something wrong in atomic_make_hmat, pls check carefull"
                    endif ! back if ( sgn /= 0 ) block

                endif ! back if ((code(delta) == 1) .and. (code(gamma) == 1)) block
            endif ! back if ((code(alpha) == 0) .and. (code(betta) == 0)) block

        enddo deltaloop ! over delta={gamma+1,norbs} loop
        enddo gammaloop ! over gamma={1,norbs-1} loop
        enddo bettaloop ! over betta={alpha+1,norbs} loop
        enddo alphaloop ! over alpha={1,norbs-1} loop

     enddo ! over jbas={1,ncfgs} loop

! dump Hamiltonian matrix to file
     call atomic_dump_hmtrx(ncfgs, Hmat)

     return
  end subroutine atomic_make_hmtrx

  subroutine atomic_make_npmat(norbs, ncfgs, invcd, npmat)
     use constants

     implicit none

! external arguments
! number of orbits
     integer, intent(in) :: norbs

! number of configurations
     integer, intent(in) :: ncfgs

! binary representation of Fock basis
     integer, intent(in) :: invcd(norbs, ncfgs)

! nn operator matrix in many body basis
     complex(dp), intent(out) :: npmat(ncfgs, ncfgs)

! local variables
! loop index over orbits
     integer :: iorb

! loop index over configurations
     integer :: ibas

! initialize nnmat to be zero
     npmat = czero

     do ibas=1,ncfgs
         do iorb=1,norbs
             if (invcd(iorb, ibas) .eq. 1) then
                 npmat(ibas, ibas) = npmat(ibas, ibas) + done
             endif
         enddo ! over iorb={1,norbs} loop
     enddo ! over ibas={1,ncfgs} loop

     return
  end subroutine atomic_make_npmat

! make single occupied number in many body space ref:  prb_86.155158(2012).added by syp
  subroutine atomic_make_sgnpmat(norbs, ncfgs, invcd, sgnpmat)
     use constants

     implicit none

! external arguments
! number of orbits
     integer, intent(in) :: norbs

! number of configurations
     integer, intent(in) :: ncfgs

! binary representation of Fock basis
     integer, intent(in) :: invcd(norbs, ncfgs)

! nn operator matrix in many body basis
     complex(dp), intent(out) :: sgnpmat(ncfgs, ncfgs)

! local variables
! loop index over orbits
     integer :: iorb

! loop index over configurations
     integer :: ibas

! initialize nnmat to be zero
     sgnpmat = czero

     do ibas=1,ncfgs
         do iorb=1,norbs-1,2
             if ((invcd(iorb, ibas) .eq. 1 .and. invcd(iorb+1, ibas) .eq. 0) &
                 .or. (invcd(iorb, ibas) .eq. 0 .and. invcd(iorb+1, ibas) .eq. 1) ) then
                 sgnpmat(ibas, ibas) = sgnpmat(ibas, ibas) + done
             endif
         enddo ! over iorb={1,norbs} loop
     enddo ! over ibas={1,ncfgs} loop

     return
  end subroutine atomic_make_sgnpmat

!>>> diagonalize atomic Hamiltonian matrix by using {jz} as good quantum number
!  subroutine atomic_diag_Hmat3(norbs, ncfgs, nstat, npmat, sgnpmat, c4mat, Amtrx, Bmtrx, Hmtrx, neigs, c4eig, Aeigs, Beigs, Heigs, Heigv)
  subroutine atomic_diag_Hmat3(norbs, ncfgs, nstat, npmat, sgnpmat, c4mat, Amtrx, Bmtrx, Hmtrx, neigs, Aeigs, Beigs, Heigs, Heigv)
     use constants

     implicit none

! number of orbits
     integer, intent(in) :: norbs

! number of configurations
     integer, intent(in) :: ncfgs

! number of configuration in each itot subspace
     integer, intent(in) :: nstat(0:norbs)

! single particle operator matrix in Fock space
     complex(dp), intent(in) :: npmat(ncfgs,ncfgs)

! single particle operator matrix in Fock space
     complex(dp), intent(in) :: sgnpmat(ncfgs,ncfgs)

! single particle operator matrix in Fock space
     complex(dp), intent(inout) :: c4mat(ncfgs,ncfgs)

! symmetry operator matrix in many body basis
     complex(dp), intent(in) :: Amtrx(ncfgs, ncfgs)

! symmetry operator matrix im many body bosis
     complex(dp), intent(in) :: Bmtrx(ncfgs, ncfgs)

! Hamiltonian operator matrix in many body basis
     complex(dp), intent(in) :: Hmtrx(ncfgs, ncfgs)

!the N good quantums number of every eigenstates     
     real(8), intent(out) :: neigs(ncfgs)

!the c4 good quantums number of every eigenstates     
!     real(8), intent(out) :: c4eig(ncfgs)

! eigenvalue of A operator
     real(8), intent(out) :: Aeigs(ncfgs)

! eigenvalue of A operator
     real(8), intent(out) :: Beigs(ncfgs)

! eigenvalue of atomic Hamiltonian
     real(dp), intent(out) :: Heigs(ncfgs)

! eigenvector of atomic hamiltonian
     complex(dp), intent(out) :: Heigv(ncfgs, ncfgs)

! local variables
! loop index over configurations
     integer :: ibas, jbas

!     real(8) :: neigs(ncfgs)

! number of dimensions in each subspace
     integer :: nbase, ndims
     integer :: is, ie, i, j

! status while allocating memory
     integer :: istat

! auxiliary integer variables and matrix
     integer :: itmp,syp_flag
     integer, allocatable :: Adegs(:)
     integer, allocatable :: Bdegs(:)
     integer, allocatable :: Hdegs(:)
     integer, allocatable :: c4degs(:)
     integer, allocatable :: sgdegs(:)

! auxiliary real(dp) matrix
     real(dp) :: jzval, jjval
     real(dp) :: Avals, Bvals, dtmp, Nvals, c4vals

! auxiliary complex(dp) matrix
     complex(dp), allocatable :: Aeigv(:, :)
     complex(dp), allocatable :: Beigv(:, :)
     complex(dp), allocatable :: c4eigv(:, :)
     complex(dp), allocatable :: sgeigv(:, :)
     complex(dp), allocatable :: sgeigs(:)

     complex(dp), allocatable :: zeigv(:, :)
     complex(dp), allocatable :: zmata(:, :)
     complex(dp), allocatable :: zmatb(:, :)
     complex(dp), allocatable :: c4mata(:, :)
     complex(dp), allocatable :: c4matb(:, :)
     complex(dp), allocatable :: tranmat(:, :)
     complex(dp), allocatable :: tranmat1(:, :)
  

! allocate memory for auxiliary arrays
     allocate(zeigv(ncfgs, ncfgs), stat=istat)
     allocate(zmata(ncfgs, ncfgs), stat=istat)
     allocate(zmatb(ncfgs, ncfgs), stat=istat)
     allocate(c4mata(ncfgs, ncfgs), stat=istat)
     allocate(c4matb(ncfgs, ncfgs), stat=istat)
     allocate(tranmat(ncfgs, ncfgs), stat=istat)
     allocate(tranmat1(ncfgs, ncfgs), stat=istat)
     allocate(Adegs(ncfgs), Aeigv(ncfgs, ncfgs), stat=istat)
     allocate(Bdegs(ncfgs), Beigv(ncfgs, ncfgs), stat=istat)
     allocate(Hdegs(ncfgs), stat=istat)
     allocate(c4degs(ncfgs), c4eigv(ncfgs, ncfgs), stat=istat)
     allocate(sgdegs(ncfgs), sgeigv(ncfgs, ncfgs), stat=istat)
     allocate(sgeigs(ncfgs), stat=istat)
     if (istat .ne. 0) stop "allocate memory error in atomic_diag_Hmat3"

! setup degeneracy structure of jzeig
     syp_flag = 1
     tranmat  = dcmplx(0.d0, 0.d0)
     tranmat1 = dcmplx(0.d0, 0.d0)
     do i = 1, ncfgs
         tranmat(i, i) = dcmplx(1.d0, 0.d0)
         tranmat1(i, i) = dcmplx(1.d0, 0.d0)
     enddo
     if(syp_flag .eq. 1)then
         call zmat_zbdev3(norbs, ncfgs, nstat,  5, 7, sgnpmat, sgeigs, sgeigv)
         call atomic_make_edeg3(norbs, ncfgs, nstat, 5, 7, sgeigs, sgdegs)
!         call zmat_zbdev3(norbs, ncfgs,nstat,  5, 7, Amtrx, Aeigs, Aeigv)
!         call atomic_make_edeg3(norbs, ncfgs, nstat, 5, 7, Aeigs, Adegs)
     elseif(ncfgs .lt. sum(Nstat)) then
         call zmat_zheev(ncfgs, ncfgs, Amtrx, Aeigs, Aeigv)
         call atomic_make_edeg0(ncfgs, Aeigs, Adegs)
     else 
         call zmat_zbdev0(norbs, ncfgs, Nstat, Amtrx, Aeigs, Aeigv)
         call atomic_make_edegn(norbs, ncfgs, Nstat, Aeigs, Adegs)
     endif ! back if (ncfgs .lt. sum(nstat)) block

! -------------------------
! transform jz operator matrix into sg operator's eigen-basis
     call atomic_dump_c4mat2(ncfgs, sgeigv,'test_sg_eigv.dat')
     call atomic_dump_isgnmat(ncfgs, sgdegs, 'test_sgnmat_degs.dat')
     call zmat_zgemm0(ncfgs, Amtrx, sgeigv, zmatb)
     call zmat_zgemm2(ncfgs, sgeigv, zmatb, zmata)
     call atomic_dump_c4mat2(ncfgs, zmata,'test_jzmat_sg.dat')


! take into account of block diagonal structure of jz operator matrix
     call zmat_zbdev1(ncfgs, sgdegs, zmata, Aeigs, Aeigv)
     call atomic_make_edeg1(ncfgs, sgdegs, Aeigs, Adegs)
     call atomic_dump_isgnmat(ncfgs, Adegs, 'test_Sz_degs.dat')
     call atomic_dump_rsgnmat(ncfgs, Aeigs, 'test_Sz_eigs.dat')
     call atomic_dump_c4mat2(ncfgs, Aeigv, "test_sz_eigv.dat")
! --------------------------

! transform jj operator matrix into sg operator's eigen-basis
     call zmat_zgemm0(ncfgs, Bmtrx, sgeigv, zmatb)
     call zmat_zgemm2(ncfgs, sgeigv, zmatb, zmata)
     call atomic_dump_c4mat2(ncfgs, zmata,'test_jjmat_sg.dat')

! transform jj operator matrix into jz operator's eigen-basis
     call zmat_zgemm0(ncfgs, zmata, Aeigv, zmatb)
     call zmat_zgemm2(ncfgs, Aeigv, zmatb, zmata)
     call atomic_dump_c4mat2(ncfgs, zmata,'test_jjmat_sz.dat')

! take into account of block diagonal structure of j2 operator matrix
     call zmat_zbdev1(ncfgs, Adegs, zmata, Beigs, Beigv)
     call atomic_make_edeg1(ncfgs, Adegs, Beigs, Bdegs)
     call atomic_dump_isgnmat(ncfgs, Bdegs, 'test_jj_degs.dat')
! -------------------------

! transform c4 operator matrix into sg operator's eigen-basis
     call zmat_zgemm0(ncfgs, c4mat, sgeigv, c4matb)
     call zmat_zgemm2(ncfgs, sgeigv, c4matb, c4mata)
     call atomic_dump_c4mat2(ncfgs, c4mata,'test_c4_sg.dat')

! transform c4 operator matrix into jz operator's eigen-basis
     call zmat_zgemm0(ncfgs, c4mata, Aeigv, c4matb)
     call zmat_zgemm2(ncfgs, Aeigv, c4matb, c4mata)
     call atomic_dump_c4mat2(ncfgs, c4mata,'test_c4_jz.dat')

! transform c4 operator matrix into j2 operator's eigen-basis
     call zmat_zgemm0(ncfgs, c4mata, Beigv, c4matb)
     call zmat_zgemm2(ncfgs, Beigv, c4matb, c4mata)
     call atomic_dump_c4mat2(ncfgs, c4mata,'test_c4_jj.dat')
! ---------------------------------


     call atomic_dump_c4mat2(ncfgs, Hmtrx,'test_ham_0.dat')
! transform Hamiltonian matrix into sg operator's eigen-basis
     call zmat_zgemm0(ncfgs, Hmtrx, sgeigv, zmatb)
     call zmat_zgemm2(ncfgs, sgeigv, zmatb, zmata)
     call atomic_dump_c4mat2(ncfgs, zmata,'test_ham_sg.dat')

! transform Hamiltonian matrix into jz operator's eigen-basis
     call zmat_zgemm0(ncfgs, zmata, Aeigv, zmatb)
     call zmat_zgemm2(ncfgs, Aeigv, zmatb, zmata)
     call atomic_dump_c4mat2(ncfgs, zmata,'test_ham_jz.dat')

! transform Hamiltonian matrix into j2 operator's eigen-basis
     call zmat_zgemm0(ncfgs, zmata, Beigv, zmatb)
     call zmat_zgemm2(ncfgs, Beigv, zmatb, zmata)
     call atomic_dump_c4mat2(ncfgs, zmata,'test_Ham_jj.dat')

! take into account of block diagonal structure of Hamitonian matrix
     call zmat_zbdev1(ncfgs, Bdegs, zmata, Heigs, zeigv)
     call atomic_make_edeg1(ncfgs, Bdegs, Heigs, Hdegs)
     call atomic_dump_isgnmat(ncfgs, Hdegs, 'test_Ham_degs.dat')
     call atomic_dump_rsgnmat(ncfgs, Heigs, 'test_Ham_eigs.dat')
     call atomic_dump_c4mat2(ncfgs, zeigv,'test_ham_vec.dat')
     
! ----------------------

! transfrom eigvectors back to orginal fock basis
     call zmat_zgemm0(ncfgs, sgeigv, Aeigv, zmata)
     call zmat_zgemm0(ncfgs, zmata,  Beigv, zmatb)
     call zmat_zgemm0(ncfgs, zmatb,  zeigv, Heigv)
     call atomic_dump_c4mat2(ncfgs,  Heigv,'test_Heigv_fock.dat')
     

! transform c4 operator into Heigv representation
     call zmat_zgemm0(ncfgs, c4mata, zeigv, c4matb)
     call zmat_zgemm2(ncfgs, zeigv, c4matb, c4mata)
     call atomic_dump_c4mat2(ncfgs, c4mata,'test_c4_En.dat')

! label every eigenstates with the single particle number
     call atomic_make_neigs(ncfgs, npmat, Heigv, neigs) 
     call atomic_dump_rsgnmat(ncfgs, neigs, 'test_npeig_ham.dat')

!---------------------------------------------------------------------!
! use selection sort to minimize swaps of eigenvectors, ref: dsteqr.f !
!---------------------------------------------------------------------!
!$   do ibas=1,ncfgs-1
!$       itmp = ibas; dtmp = eigval(itmp)

!$       do jbas=ibas+1,ncfgs
!$           if (eigval(jbas) .lt. dtmp) then
!$               itmp = jbas; dtmp = eigval(jbas)
!$           endif
!$       enddo ! over jbas={ibas+1,ncfgs} loop

!$       if (itmp .ne. ibas) then
!$           eigval(itmp) = eigval(ibas); eigval(ibas) = dtmp
!$           call zswap(ncfgs, eigvec(1, ibas), 1, eigvec(1,itmp), 1)
!$       endif
!$   enddo ! over ibas={2,ncfgs} loop

!---------------------------------------------------------------------!
! sort energy ascending order, jzeig in descending order
!---------------------------------------------------------------------!
     do ibas=1,ncfgs-1
         itmp  = ibas        ; Bvals = Beigs(itmp)
         dtmp  = Heigs(itmp) ; Avals = Aeigs(itmp)
         Nvals = neigs(itmp)
         tranmat1 = dcmplx(0.d0, 0.d0)
         do i = 1, ncfgs
             tranmat1(i, i) = dcmplx(1.d0, 0.d0)
         enddo
         do jbas=ibas+1,ncfgs
             if(Heigs(jbas) + 1.0d-6 .lt. dtmp)then
                 itmp  = jbas        ; Bvals = Beigs(jbas)
                 Nvals = neigs(jbas) ; Avals = Aeigs(jbas)
                 dtmp  = Heigs(jbas)
             else if (abs(Heigs(jbas)-dtmp) .lt. 1.0d-6)then
                 if (neigs(jbas) + 1.0d-12 .lt. Nvals) then
                     itmp = jbas         ; Bvals = Beigs(jbas)
                     Nvals = neigs(jbas) ; Avals = Aeigs(jbas)
                 else if (abs(neigs(jbas)-Nvals) .lt. 1.0d-12) then
                     if (Beigs(jbas) + 1.0d-4 .lt. Bvals) then
                         itmp = jbas
                         Bvals = Beigs(jbas); Avals = Aeigs(jbas)
                     else if (abs(Beigs(jbas)-Bvals) .lt. 1.0d-4) then
                         if (Aeigs(jbas) .gt. Avals) then
                             itmp = jbas; Avals = Aeigs(jbas)
                         endif
                     endif
                 endif
             endif
         enddo ! over jbas={ibas+1,ncfgs} loop

         if (itmp .ne. ibas) then
             Heigs(itmp) = Heigs(ibas); Heigs(ibas) = dtmp
             Aeigs(itmp) = Aeigs(ibas); Aeigs(ibas) = Avals
             Beigs(itmp) = Beigs(ibas); Beigs(ibas) = Bvals
             neigs(itmp) = neigs(ibas); neigs(ibas) = Nvals
             call zswap(ncfgs, Heigv(1, ibas), 1, Heigv(1,itmp), 1)
             tranmat1(itmp, itmp) = dcmplx(0.d0, 0.d0)
             tranmat1(itmp, ibas) = dcmplx(1.d0, 0.d0)
             tranmat1(ibas, itmp) = dcmplx(1.d0, 0.d0)
             tranmat1(ibas, ibas) = dcmplx(0.d0, 0.d0)
             tranmat = matmul(tranmat, tranmat1)
         endif

     enddo ! over ibas={2,ncfgs} loop
     ! transform c4 operator into new energy eigenstates
     call zmat_zgemm0(ncfgs, c4mata, tranmat, c4matb)
     call zmat_zgemm2(ncfgs, tranmat, c4matb, c4mat)
     call atomic_dump_c4mat2(ncfgs, c4mat,'test_c4_newEn.dat')
!---------------------------------------------------------------------!
!     do ibas=1,ncfgs-1
!         itmp = ibas        ; Bvals = Beigs(itmp)
!         dtmp = Heigs(itmp) ; Avals = Aeigs(itmp)
!
!         do jbas=ibas+1,ncfgs
!             if (Heigs(jbas) + 1.0d-12 .lt. dtmp) then
!                 itmp = jbas       ; Bvals = Beigs(jbas)
!                 dtmp = Heigs(jbas); Avals = Aeigs(jbas)
!             else if (abs(Heigs(jbas)-dtmp) .lt. 1.0d-12) then
!                 if (Beigs(jbas) + 1.0d-4 .lt. Bvals) then
!                     itmp = jbas
!                     Bvals = Beigs(jbas); Avals = Aeigs(jbas)
!                 else if (abs(Beigs(jbas)-Bvals) .lt. 1.0d-4) then
!                     if (Aeigs(jbas) .gt. Avals) then
!                         itmp = jbas; Avals = Aeigs(jbas)
!                     endif
!                 endif
!             endif
!         enddo ! over jbas={ibas+1,ncfgs} loop
!
!         if (itmp .ne. ibas) then
!             Heigs(itmp) = Heigs(ibas); Heigs(ibas) = dtmp
!             Aeigs(itmp) = Aeigs(ibas); Aeigs(ibas) = Avals
!             Beigs(itmp) = Beigs(ibas); Beigs(ibas) = Bvals
!             call zswap(ncfgs, Heigv(1, ibas), 1, Heigv(1,itmp), 1)
!         endif
!
!     enddo ! over ibas={2,ncfgs} loop

     return
  end subroutine atomic_diag_Hmat3


!>>> diagonalize atomic Hamiltonian matrix by using {N, N_single,jz} as good quantum number in every irrep basis block
!  subroutine atomic_diag_Hmat3(norbs, ncfgs, nstat, npmat, sgnpmat, c4mat, Amtrx, Bmtrx, Hmtrx, neigs, c4eig, Aeigs, Beigs, Heigs, Heigv)
  subroutine atomic_diag_Hmat4(norbs, ncfgs, nstat, c4virrepdegs, Umtrx, npmat, sgnpmat,  &
          Amtrx, Bmtrx, Hmtrx, neigs, sgeigs, Aeigs, Beigs, Heigs, Heigv, Hdegs)
     use constants

     implicit none

! number of orbits
     integer, intent(in) :: norbs

! number of configurations
     integer, intent(in) :: ncfgs

! number of configuration in each itot subspace
     integer, intent(in) :: nstat(0:norbs)

! the dimensionility of every block of irrep basis of C4v
     integer, intent(in) :: c4virrepdegs(ncfgs)

! single particle operator matrix in Fock space
     complex(dp), intent(in) :: Umtrx(ncfgs,ncfgs)

! single particle operator matrix in Fock space
     complex(dp), intent(in) :: npmat(ncfgs,ncfgs)

! single particle operator matrix in Fock space
     complex(dp), intent(in) :: sgnpmat(ncfgs,ncfgs)

! symmetry operator matrix in many body basis
     complex(dp), intent(in) :: Amtrx(ncfgs, ncfgs)

! symmetry operator matrix im many body bosis
     complex(dp), intent(in) :: Bmtrx(ncfgs, ncfgs)

! Hamiltonian operator matrix in many body basis
     complex(dp), intent(in) :: Hmtrx(ncfgs, ncfgs)

!the N good quantums number of every eigenstates     
     real(8), intent(out) :: neigs(ncfgs)

!the N_{sg} good quantums number of every eigenstates     
     real(8), intent(out) :: sgeigs(ncfgs)

!the c4 good quantums number of every eigenstates     
!     real(8), intent(out) :: c4eig(ncfgs)

! eigenvalue of A operator
     real(8), intent(out) :: Aeigs(ncfgs)

! eigenvalue of A operator
     real(8), intent(out) :: Beigs(ncfgs)

! eigenvalue of atomic Hamiltonian
     real(dp), intent(out) :: Heigs(ncfgs)

! eigenvector of atomic hamiltonian
     complex(dp), intent(out) :: Heigv(ncfgs, ncfgs)

! degeneracy of atomic eigenvalues
     integer,     intent(out) :: Hdegs(ncfgs)

! local variables
! loop index over configurations
     integer :: ibas, jbas

!     real(8) :: neigs(ncfgs)

! number of dimensions in each subspace
     integer :: nbase, ndims, num_vpmstot
     integer :: is, ie, i, j

! status while allocating memory
     integer :: istat

! auxiliary integer variables and matrix
     integer :: itmp,syp_flag
     integer, allocatable :: Adegs(:)
     integer, allocatable :: Bdegs(:)
!     integer, allocatable :: Hdegs(:)
     integer, allocatable :: c4degs(:)
     integer, allocatable :: sgdegs(:)
     integer, allocatable :: Ndegs(:)
     integer, allocatable :: list_deg(:)

! auxiliary real(dp) matrix
     real(dp) :: jzval, jjval
     real(dp) :: Avals, Bvals, dtmp, Nvals, c4vals

! auxiliary complex(dp) matrix
     complex(dp), allocatable :: Aeigv(:, :)
     complex(dp), allocatable :: Beigv(:, :)
     complex(dp), allocatable :: c4eigv(:, :)
     complex(dp), allocatable :: sgeigv(:, :)
!     real(dp), allocatable :: sgeigs(:)
     complex(dp), allocatable :: Neigv(:, :)
     real(dp), allocatable :: Neigs1(:)

     complex(dp), allocatable :: zeigv(:, :)
     complex(dp), allocatable :: zmata(:, :)
     complex(dp), allocatable :: zmatb(:, :)
     complex(dp), allocatable :: c4mata(:, :)
     complex(dp), allocatable :: c4matb(:, :)
     complex(dp), allocatable :: tranmat(:, :)
     complex(dp), allocatable :: tranmat1(:, :)
  

! allocate memory for auxiliary arrays
     allocate(zeigv(ncfgs, ncfgs), stat=istat)
     allocate(zmata(ncfgs, ncfgs), stat=istat)
     allocate(zmatb(ncfgs, ncfgs), stat=istat)
     allocate(c4mata(ncfgs, ncfgs), stat=istat)
     allocate(c4matb(ncfgs, ncfgs), stat=istat)
     allocate(tranmat(ncfgs, ncfgs), stat=istat)
     allocate(tranmat1(ncfgs, ncfgs), stat=istat)
     allocate(Adegs(ncfgs), Aeigv(ncfgs, ncfgs), stat=istat)
     allocate(Bdegs(ncfgs), Beigv(ncfgs, ncfgs), stat=istat)
!     allocate(Hdegs(ncfgs), stat=istat)
     allocate(c4degs(ncfgs), c4eigv(ncfgs, ncfgs), stat=istat)
     allocate(sgdegs(ncfgs), sgeigv(ncfgs, ncfgs), stat=istat)
!     allocate(sgeigs(ncfgs), stat=istat)
     allocate(Ndegs(ncfgs), stat=istat)
     allocate(Neigv(ncfgs, ncfgs), stat=istat)
     allocate(Neigs1(ncfgs), stat=istat)
     allocate(list_deg(ncfgs), stat=istat)
     if (istat .ne. 0) stop "allocate memory error in atomic_diag_Hmat3"

! setup degeneracy structure of jzeig
     syp_flag = 1
     tranmat  = dcmplx(0.d0, 0.d0)
     tranmat1 = dcmplx(0.d0, 0.d0)
     do i = 1, ncfgs
         tranmat(i, i) = dcmplx(1.d0, 0.d0)
         tranmat1(i, i) = dcmplx(1.d0, 0.d0)
     enddo
!     if(syp_flag .eq. 1)then
!         call zmat_zbdev3(norbs, ncfgs, nstat,  5, 7, sgnpmat, sgeigs, sgeigv)
!         call atomic_make_edeg3(norbs, ncfgs, nstat, 5, 7, sgeigs, sgdegs)
!!         call zmat_zbdev3(norbs, ncfgs,nstat,  5, 7, Amtrx, Aeigs, Aeigv)
!!         call atomic_make_edeg3(norbs, ncfgs, nstat, 5, 7, Aeigs, Adegs)
!     elseif(ncfgs .lt. sum(Nstat)) then
!         call zmat_zheev(ncfgs, ncfgs, Amtrx, Aeigs, Aeigv)
!         call atomic_make_edeg0(ncfgs, Aeigs, Adegs)
!     else 
!         call zmat_zbdev0(norbs, ncfgs, Nstat, Amtrx, Aeigs, Aeigv)
!         call atomic_make_edegn(norbs, ncfgs, Nstat, Aeigs, Adegs)
!     endif ! back if (ncfgs .lt. sum(nstat)) block

! -------------------------
     call count_vpm(ncfgs, c4virrepdegs, num_vpmstot, list_deg)
     call atomic_dump_deglist(ncfgs, num_vpmstot, list_deg, "c4v_deglist.dat")

! -------------------------
! take into account of block diagonal structure of npmat

     call zmat_zbdev1(ncfgs, c4virrepdegs, npmat, neigs1, neigv)
     call atomic_make_edeg1(ncfgs, c4virrepdegs, neigs1, ndegs)
     call atomic_dump_c4mat2(ncfgs, neigv,'test_npmat_eigv.c4v.dat')
     call atomic_dump_isgnmat(ncfgs, ndegs, 'test_npmat_degs.dat')
     call count_vpm(ncfgs, ndegs, num_vpmstot, list_deg)
     write(20,*)'npmat_deglist', list_deg
     call atomic_dump_deglist(ncfgs, num_vpmstot, list_deg, "npmat_deglist.dat")

! ------------------------


! transform "sgnmat" operator matrix into particle number operator's eigen-basis
     call zmat_zgemm0(ncfgs, sgnpmat, neigv, zmatb)
     call zmat_zgemm2(ncfgs, neigv, zmatb, zmata)
     call atomic_dump_c4mat2(ncfgs, zmata,'test_sgnpmat.n.dat')


! take into account of block diagonal structure of "sgnmat" operator matrix
     call zmat_zbdev1(ncfgs, ndegs, zmata, sgeigs, sgeigv)
     call atomic_make_edeg1(ncfgs, ndegs, sgeigs, sgdegs)
     call atomic_dump_isgnmat(ncfgs, sgdegs, 'test_sg_degs.dat')
     call atomic_dump_rsgnmat(ncfgs, sgeigs, 'test_sg_eigs.dat')
     call atomic_dump_c4mat2(ncfgs, sgeigv, "test_sg_eigv.n.dat")
     call count_vpm(ncfgs, sgdegs, num_vpmstot, list_deg)
     call atomic_dump_deglist(ncfgs, num_vpmstot, list_deg, "sgnpmat_deglist.dat")
! --------------------------

! transform ss operator matrix into particle number operator's eigen-basis
     call zmat_zgemm0(ncfgs, Bmtrx, neigv, zmatb)
     call zmat_zgemm2(ncfgs, neigv, zmatb, zmata)
     call atomic_dump_c4mat2(ncfgs, zmata,'test_ssmat.n.dat')


! transform jz operator matrix into sg operator's eigen-basis
     call zmat_zgemm0(ncfgs, zmata, sgeigv, zmatb)
     call zmat_zgemm2(ncfgs, sgeigv, zmatb, zmata)
     call atomic_dump_c4mat2(ncfgs, zmata,'test_ssmat.sg.dat')


! take into account of block diagonal structure of ss operator matrix
     call zmat_zbdev1(ncfgs, sgdegs, zmata, Beigs, Beigv)
     call atomic_make_edeg1(ncfgs, sgdegs, Beigs, Bdegs)
     call atomic_dump_isgnmat(ncfgs, Bdegs, 'test_ss_degs.dat')
     call atomic_dump_rsgnmat(ncfgs, Beigs, 'test_ss_eigs.dat')
     call atomic_dump_c4mat2(ncfgs, Beigv, "test_ss_eigv.sg.dat")
     call count_vpm(ncfgs, Bdegs, num_vpmstot, list_deg)
     call atomic_dump_deglist(ncfgs, num_vpmstot, list_deg, "ss_deglist.dat")
! --------------------------

! transform sz operator matrix into particle number operator's eigen-basis
     call zmat_zgemm0(ncfgs, Amtrx, neigv, zmatb)
     call zmat_zgemm2(ncfgs, neigv, zmatb, zmata)
     call atomic_dump_c4mat2(ncfgs, zmata,'test_sz.n.dat')

! transform sz operator matrix into sg operator's eigen-basis
     call zmat_zgemm0(ncfgs, zmata, sgeigv, zmatb)
     call zmat_zgemm2(ncfgs, sgeigv, zmatb, zmata)
     call atomic_dump_c4mat2(ncfgs, zmata,'test_sz.sg.dat')

! transform jj operator matrix into jz operator's eigen-basis
     call zmat_zgemm0(ncfgs, zmata, Beigv, zmatb)
     call zmat_zgemm2(ncfgs, Beigv, zmatb, zmata)
     call atomic_dump_c4mat2(ncfgs, zmata,'test_sz.ss.dat')

! take into account of block diagonal structure of j2 operator matrix
     call zmat_zbdev1(ncfgs, Bdegs, zmata, Aeigs, Aeigv)
     call atomic_make_edeg1(ncfgs, Bdegs, Aeigs, Adegs)
     call atomic_dump_isgnmat(ncfgs, Adegs, 'test_sz_degs.dat')
     call count_vpm(ncfgs, Adegs, num_vpmstot, list_deg)
     call atomic_dump_deglist(ncfgs, num_vpmstot, list_deg, "sz_deglist.dat")
! -------------------------

!! transform c4 operator matrix into sg operator's eigen-basis
!     call zmat_zgemm0(ncfgs, c4mat, sgeigv, c4matb)
!     call zmat_zgemm2(ncfgs, sgeigv, c4matb, c4mata)
!     call atomic_dump_c4mat2(ncfgs, c4mata,'test_c4_sg.dat')

! transform c4 operator matrix into jz operator's eigen-basis
!     call zmat_zgemm0(ncfgs, c4mata, Aeigv, c4matb)
!     call zmat_zgemm2(ncfgs, Aeigv, c4matb, c4mata)
!     call atomic_dump_c4mat2(ncfgs, c4mata,'test_c4_jz.dat')

! transform c4 operator matrix into j2 operator's eigen-basis
!     call zmat_zgemm0(ncfgs, c4mata, Beigv, c4matb)
!     call zmat_zgemm2(ncfgs, Beigv, c4matb, c4mata)
!     call atomic_dump_c4mat2(ncfgs, c4mata,'test_c4_jj.dat')
! ---------------------------------

! ---------------------------------
! transform Hamiltonian matrix into particle number operator's eigen-basis
     call zmat_zgemm0(ncfgs, Hmtrx, neigv, zmatb)
     call zmat_zgemm2(ncfgs, neigv, zmatb, zmata)
     call atomic_dump_c4mat2(ncfgs, zmata,'test_ham.n.dat')

! transform Hamiltonian matrix into sg operator's eigen-basis
     call zmat_zgemm0(ncfgs, zmata, sgeigv, zmatb)
     call zmat_zgemm2(ncfgs, sgeigv, zmatb, zmata)
     call atomic_dump_c4mat2(ncfgs, zmata,'test_ham.sg.dat')

! transform Hamiltonian matrix into ss operator's eigen-basis
     call zmat_zgemm0(ncfgs, zmata, Beigv, zmatb)
     call zmat_zgemm2(ncfgs, Beigv, zmatb, zmata)
     call atomic_dump_c4mat2(ncfgs, zmata,'test_ham.ss.dat')

! transform Hamiltonian matrix into sz operator's eigen-basis
     call zmat_zgemm0(ncfgs, zmata, Aeigv, zmatb)
     call zmat_zgemm2(ncfgs, Aeigv, zmatb, zmata)
     call atomic_dump_c4mat2(ncfgs, zmata,'test_Ham.sz.dat')

! take into account of block diagonal structure of Hamitonian matrix
     call zmat_zbdev1(ncfgs, Adegs, zmata, Heigs, zeigv)
     call atomic_make_edeg1(ncfgs, Adegs, Heigs, Hdegs)
     call atomic_dump_isgnmat(ncfgs, Hdegs, 'test_Ham_degs.dat')
     call atomic_dump_rsgnmat(ncfgs, Heigs, 'test_Ham_eigs.dat')
     call atomic_dump_c4mat2(ncfgs, zeigv,'test_ham_vec.sz.dat')
     call count_vpm(ncfgs, Hdegs, num_vpmstot, list_deg)
     call atomic_dump_deglist(ncfgs, num_vpmstot, list_deg, "Ham_deglist.dat")
     
! ----------------------

! transfrom eigvectors back to orginal fock basis
     call zmat_zgemm0(ncfgs, Umtrx,  neigv,  zmata)
     call zmat_zgemm0(ncfgs, zmata,  sgeigv, zmatb)
     call zmat_zgemm0(ncfgs, zmatb, Beigv,  zmata)
     call zmat_zgemm0(ncfgs, zmata,  Aeigv,  zmatb)
     call zmat_zgemm0(ncfgs, zmatb,  zeigv,  Heigv)
     call atomic_dump_c4mat2(ncfgs,  Heigv,'test_Heigv.fock.dat')
     
!! transform c4 operator into Heigv representation
!     call zmat_zgemm0(ncfgs, c4mata, zeigv, c4matb)
!     call zmat_zgemm2(ncfgs, zeigv, c4matb, c4mata)
!     call atomic_dump_c4mat2(ncfgs, c4mata,'test_c4_En.dat')

! label every eigenstates with the single particle number
!     call atomic_make_neigs(ncfgs, npmat, Heigv, neigs) 
!     call atomic_dump_rsgnmat(ncfgs, neigs, 'test_npeig_ham.dat')
     call atomic_dump_rsgnmat(ncfgs, neigs1, 'test_npeig1_ham.dat')
     neigs = neigs1
!---------------------------------------------------------------------!
! use selection sort to minimize swaps of eigenvectors, ref: dsteqr.f !
!---------------------------------------------------------------------!
!$   do ibas=1,ncfgs-1
!$       itmp = ibas; dtmp = eigval(itmp)

!$       do jbas=ibas+1,ncfgs
!$           if (eigval(jbas) .lt. dtmp) then
!$               itmp = jbas; dtmp = eigval(jbas)
!$           endif
!$       enddo ! over jbas={ibas+1,ncfgs} loop

!$       if (itmp .ne. ibas) then
!$           eigval(itmp) = eigval(ibas); eigval(ibas) = dtmp
!$           call zswap(ncfgs, eigvec(1, ibas), 1, eigvec(1,itmp), 1)
!$       endif
!$   enddo ! over ibas={2,ncfgs} loop

!---------------------------------------------------------------------!
! sort energy ascending order, jzeig in descending order
!---------------------------------------------------------------------!
!     do ibas=1,ncfgs-1
!         itmp  = ibas        ; Bvals = Beigs(itmp)
!         dtmp  = Heigs(itmp) ; Avals = Aeigs(itmp)
!         Nvals = neigs(itmp)
!         tranmat1 = dcmplx(0.d0, 0.d0)
!         do i = 1, ncfgs
!             tranmat1(i, i) = dcmplx(1.d0, 0.d0)
!         enddo
!         do jbas=ibas+1,ncfgs
!             if(Heigs(jbas) + 1.0d-6 .lt. dtmp)then
!                 itmp  = jbas        ; Bvals = Beigs(jbas)
!                 Nvals = neigs(jbas) ; Avals = Aeigs(jbas)
!                 dtmp  = Heigs(jbas)
!             else if (abs(Heigs(jbas)-dtmp) .lt. 1.0d-6)then
!                 if (neigs(jbas) + 1.0d-12 .lt. Nvals) then
!                     itmp = jbas         ; Bvals = Beigs(jbas)
!                     Nvals = neigs(jbas) ; Avals = Aeigs(jbas)
!                 else if (abs(neigs(jbas)-Nvals) .lt. 1.0d-12) then
!                     if (Beigs(jbas) + 1.0d-4 .lt. Bvals) then
!                         itmp = jbas
!                         Bvals = Beigs(jbas); Avals = Aeigs(jbas)
!                     else if (abs(Beigs(jbas)-Bvals) .lt. 1.0d-4) then
!                         if (Aeigs(jbas) .gt. Avals) then
!                             itmp = jbas; Avals = Aeigs(jbas)
!                         endif
!                     endif
!                 endif
!             endif
!         enddo ! over jbas={ibas+1,ncfgs} loop
!
!         if (itmp .ne. ibas) then
!             Heigs(itmp) = Heigs(ibas); Heigs(ibas) = dtmp
!             Aeigs(itmp) = Aeigs(ibas); Aeigs(ibas) = Avals
!             Beigs(itmp) = Beigs(ibas); Beigs(ibas) = Bvals
!             neigs(itmp) = neigs(ibas); neigs(ibas) = Nvals
!             call zswap(ncfgs, Heigv(1, ibas), 1, Heigv(1,itmp), 1)
!             tranmat1(itmp, itmp) = dcmplx(0.d0, 0.d0)
!             tranmat1(itmp, ibas) = dcmplx(1.d0, 0.d0)
!             tranmat1(ibas, itmp) = dcmplx(1.d0, 0.d0)
!             tranmat1(ibas, ibas) = dcmplx(0.d0, 0.d0)
!             tranmat = matmul(tranmat, tranmat1)
!         endif
!
!     enddo ! over ibas={2,ncfgs} loop
     ! transform c4 operator into new energy eigenstates
!     call zmat_zgemm0(ncfgs, c4mata, tranmat, c4matb)
!     call zmat_zgemm2(ncfgs, tranmat, c4matb, c4mat)
!     call atomic_dump_c4mat2(ncfgs, c4mat,'test_c4_newEn.dat')
!---------------------------------------------------------------------!
!     do ibas=1,ncfgs-1
!         itmp = ibas        ; Bvals = Beigs(itmp)
!         dtmp = Heigs(itmp) ; Avals = Aeigs(itmp)
!
!         do jbas=ibas+1,ncfgs
!             if (Heigs(jbas) + 1.0d-12 .lt. dtmp) then
!                 itmp = jbas       ; Bvals = Beigs(jbas)
!                 dtmp = Heigs(jbas); Avals = Aeigs(jbas)
!             else if (abs(Heigs(jbas)-dtmp) .lt. 1.0d-12) then
!                 if (Beigs(jbas) + 1.0d-4 .lt. Bvals) then
!                     itmp = jbas
!                     Bvals = Beigs(jbas); Avals = Aeigs(jbas)
!                 else if (abs(Beigs(jbas)-Bvals) .lt. 1.0d-4) then
!                     if (Aeigs(jbas) .gt. Avals) then
!                         itmp = jbas; Avals = Aeigs(jbas)
!                     endif
!                 endif
!             endif
!         enddo ! over jbas={ibas+1,ncfgs} loop
!
!         if (itmp .ne. ibas) then
!             Heigs(itmp) = Heigs(ibas); Heigs(ibas) = dtmp
!             Aeigs(itmp) = Aeigs(ibas); Aeigs(ibas) = Avals
!             Beigs(itmp) = Beigs(ibas); Beigs(ibas) = Bvals
!             call zswap(ncfgs, Heigv(1, ibas), 1, Heigv(1,itmp), 1)
!         endif
!
!     enddo ! over ibas={2,ncfgs} loop

     return
  end subroutine atomic_diag_Hmat4



! get degeneracy of eigs in each subspace of N operator
  subroutine atomic_make_edegN(norbs, ncfgs, Nstat, eigs, edeg)
     use constants

     implicit none

! number of orbits
     integer, intent(in) :: norbs

! number of configurations
     integer, intent(in) :: ncfgs

! number of configuration in each itot subspace
     integer, intent(in) :: nstat(0:norbs)

! eigenvalues of atomic Hamiltonian
     real(dp), intent(in) :: eigs(ncfgs)

! degeneracy structure of energy spectrum
     integer, intent(out) :: edeg(ncfgs)

! local variables
! loop index over configurations
     integer :: ibas
     integer :: itot

     integer :: is, ie
     integer :: nbase

! auxiliary real(dp) variables
     real(dp) :: eval

! initialize some variables
     edeg = 0; nbase = 0

     do itot=0,norbs
         eval = 1000.D0
         is = 1; ie = nstat(itot)

         do ibas=nbase+ie,nbase+is,-1
             if (abs(eigs(ibas)-eval) .lt. 1.0D-7) then
                 edeg(ibas) = edeg(ibas+1) + 1
             else
                 edeg(ibas) = 1; eval = eigs(ibas)
             endif
         enddo ! over ibas={nbase+ie,nbase+is,-1} loop

         nbase = nbase + nstat(itot)
     enddo ! over ibas={ncfgs,1,-1} loop

     return
  end subroutine atomic_make_edegN


  ! get degeneracy of eigs in selective subspace of N operator
  subroutine atomic_make_edeg3(norbs, ncfgs, Nstat, ntiny, nlarge, eigs, edeg)
     use constants

     implicit none

! number of orbits
     integer, intent(in) :: norbs

! number of configurations
     integer, intent(in) :: ncfgs

! number of configuration in each itot subspace
     integer, intent(in) :: nstat(0:norbs)

! low boundry of selected N subspace
     integer, intent(in) :: ntiny

! top boundry of selected N subspace
     integer, intent(in) :: nlarge

! eigenvalues of atomic Hamiltonian
     real(dp), intent(in) :: eigs(ncfgs)

! degeneracy structure of energy spectrum
     integer, intent(out) :: edeg(ncfgs)

! local variables
! loop index over configurations
     integer :: ibas
     integer :: itot

     integer :: is, ie
     integer :: nbase

! auxiliary real(dp) variables
     real(dp) :: eval

! initialize some variables
     edeg = 0; nbase = 0

     do itot=ntiny,nlarge
         eval = 1000.D0
         is = 1; ie = nstat(itot)

         do ibas=nbase+ie,nbase+is,-1
             if (abs(eigs(ibas)-eval) .lt. 1.0D-7) then
                 edeg(ibas) = edeg(ibas+1) + 1
             else
                 edeg(ibas) = 1; eval = eigs(ibas)
             endif
         enddo ! over ibas={nbase+ie,nbase+is,-1} loop

         nbase = nbase + nstat(itot)
     enddo ! over ibas={ncfgs,1,-1} loop

     return
  end subroutine atomic_make_edeg3


  ! get degeneracy of one of A's block representation  
  subroutine atomic_make_edeg0(ncfgs, Aeigs, Adegs)

     implicit none

! number of configurations
     integer, intent(in) :: ncfgs

! eigenvalues of atomic Hamiltonian
     real(8), intent(in) :: Aeigs(ncfgs)

! degeneracy structure of energy spectrum
     integer, intent(out) :: Adegs(ncfgs)

! local variables
! loop index over configurations
     integer :: ibas

! auxiliary real(dp) variables
     real(8) :: eval

! initialize some variables
     Adegs = 0; eval = 1.0D+7

! setup degeneracy structure of Aeigs
     do ibas=ncfgs,1,-1
         if (abs(Aeigs(ibas)-eval) .lt. 1.0D-7) then
             Adegs(ibas) = Adegs(ibas+1) + 1
         else
             Adegs(ibas) = 1; eval = Aeigs(ibas)
         endif
     enddo ! over ibas={ncfgs,1,-1} loop

     return
  end subroutine atomic_make_edeg0

! get degeneracy of B operator in A's subspace  
  subroutine atomic_make_edeg1(ncfgs, Adegs, Beigs, Bdegs)

     implicit none

! number of configurations
     integer, intent(in) :: ncfgs

     integer, intent(in) :: Adegs(ncfgs)

! eigenvalues of atomic Hamiltonian
     real(8), intent(in) :: Beigs(ncfgs)

! degeneracy structure of energy spectrum
     integer, intent(out) :: Bdegs(ncfgs)

! local variables
! loop index over configurations
     integer :: ibas

! auxiliary real(dp) variables
     real(8) :: eval

     integer :: is, ie
     integer :: nbase, ndims

! initialize some variables
     Bdegs = 0; nbase = 0

! setup degeneracy structure of Aeigs
     do while (nbase .lt.ncfgs)
         is = nbase + 1; ndims = Adegs(nbase+1)
         ie = nbase + ndims; eval = 1000000.0D0

         do ibas=ie,is,-1
             if (abs(Beigs(ibas)-eval) .lt. 1.0D-7) then
                 Bdegs(ibas) = Bdegs(ibas+1) + 1
             else
                 Bdegs(ibas) = 1; eval = Beigs(ibas)
             endif
         enddo ! over ibas={ie,is,-1} loop

         nbase = nbase + ndims
     enddo ! over ibas={ncfgs,1,-1} loop

     return
  end subroutine atomic_make_edeg1

!>>> diagonalize block-diagonal complex matrix
  subroutine zmat_zbdev1(ncfgs, ndegs, zmata, deigs, zeigv)

     implicit none

! number of configurations
     integer, intent(in) :: ncfgs

! block diagonal structure of 
     integer, intent(in) :: ndegs(ncfgs)

! input matrix for diagonalization
     complex(8), intent(in) :: zmata(ncfgs, ncfgs)

! eigenvalues of the input operator matrix
     real(8), intent(out) :: deigs(ncfgs)

! eigenvectors of the input operator matrix
     complex(8), intent(out) :: zeigv(ncfgs, ncfgs)

! offset of current deger
     integer :: nbase
     integer :: ndims
     integer :: is, ie

     deigs = 0.0D0; zeigv = dcmplx(0.0D0, 0.0D0)

     nbase = 0; ndims = ndegs(1)
     do while (nbase .lt. ncfgs)
         is = nbase + 1; ie = nbase + ndims
         call zmat_zheev(ndims, ndims, zmata(is:ie, is:ie), deigs(is:ie), zeigv(is:ie, is:ie))
         nbase = nbase + ndims; ndims = ndegs(nbase+1)
     enddo ! over while (nbase .lt. ncfgs) loop

     return
  end subroutine zmat_zbdev1

!>>> diagonalize block-diagonal complex matrix
  subroutine zmat_zbdev0(norbs, ncfgs,  nstat, zmata, deigs, zeigv)

     implicit none

! number of orbits
     integer, intent(in) :: norbs

! number of configurations
     integer, intent(in) :: ncfgs

! block diagonal structure of 
     integer, intent(in) :: nstat(0:norbs)

! input matrix for diagonalization
     complex(8), intent(in) :: zmata(ncfgs, ncfgs)

! eigenvalues of the input operator matrix
     real(8), intent(out) :: deigs(ncfgs)

! eigenvectors of the input operator matrix
     complex(8), intent(out) :: zeigv(ncfgs, ncfgs)

! offset of current deger
     integer :: nbase
     integer :: ndims
     integer :: is, ie
     integer :: itot

     deigs = 0.0D0; zeigv = dcmplx(0.0D0, 0.0D0)

     nbase = 0
     do itot=0,norbs
         ndims = nstat(itot); is = nbase + 1; ie = nbase + ndims; nbase = nbase + ndims 
         call zmat_zheev(ndims, ndims, zmata(is:ie, is:ie), deigs(is:ie), zeigv(is:ie, is:ie))
     enddo ! over while (nbase .lt. ncfgs) loop

     return
  end subroutine zmat_zbdev0


!>>> diagonalize block-diagonal complex matrix between 0 < n_1< n_2< ntots
  subroutine zmat_zbdev3(norbs, ncfgs, nstat, ntiny, nlarge, zmata, deigs, zeigv)

     implicit none

! number of orbits
     integer, intent(in) :: norbs

! number of configurations
     integer, intent(in) :: ncfgs

! block diagonal structure of 
     integer, intent(in) :: nstat(0:norbs)

! the low boundry of N subspace
     integer, intent(in) :: ntiny

! the top boundry of N subspace
     integer, intent(in) :: nlarge

! input matrix for diagonalization
     complex(8), intent(in) :: zmata(ncfgs, ncfgs)

! eigenvalues of the input operator matrix
     real(8), intent(out) :: deigs(ncfgs)

! eigenvectors of the input operator matrix
     complex(8), intent(out) :: zeigv(ncfgs, ncfgs)

! offset of current deger
     integer :: nbase
     integer :: ndims
     integer :: is, ie
     integer :: itot

     deigs = 0.0D0; zeigv = dcmplx(0.0D0, 0.0D0)

     nbase = 0
     do itot=ntiny, nlarge
         ndims = nstat(itot); is = nbase + 1; ie = nbase + ndims; nbase = nbase + ndims 
         call zmat_zheev(ndims, ndims, zmata(is:ie, is:ie), deigs(is:ie), zeigv(is:ie, is:ie))
     enddo ! over while (nbase .lt. ncfgs) loop

     return
  end subroutine zmat_zbdev3


  subroutine atomic_make_neigs(ncfgs, npmat, Heigv, neigs)
     implicit none
     integer, intent(in) :: ncfgs
     complex(8), intent(in) :: npmat(ncfgs, ncfgs)
     complex(8), intent(in) :: Heigv(ncfgs, ncfgs)
     real(8), intent(out) :: neigs(ncfgs)
     integer :: ibas

     complex(8) :: zmata(ncfgs, ncfgs), zmatb(ncfgs, ncfgs)

     call zmat_zgemm0(ncfgs, npmat, Heigv, zmatb)
     call zmat_zgemm2(ncfgs, Heigv, zmatb, zmata)

     do ibas=1,ncfgs
         neigs(ibas) = real(zmata(ibas, ibas))
     enddo

     return
  end subroutine atomic_make_neigs
