!=========================================================================!
! project : rambutan
! program : atomic_main
! history : 06/25/2011
! authors : xidai and duliang (email:duleung@gmail.com)
! purpose : solve atomic Hamiltonian with just spin-orbit coupling
!           for purpose of GVA, we setup single-particle basis in (j2, jz)
! comment : These programs are distributed in the hope that they will be 
!           useful, but WITHOUT ANY WARRANTY; without even the implied 
!           warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
!=========================================================================!
  program atomic_main
     use constants
     use control
     use context

     implicit none

! print the running header for atomic problems
     call atomic_print_header()

! setup important parameters for atomic problems
     call atomic_config()
     call atomic_dump_Upara(Uc, Uv, Jz, Js, Jp, lamb)

! print the running summary for atomic problems
     call atomic_print_summary()

! allocate memory for global arrays and initialize them
     call atomic_setup_array()

! setup key matrix in atomic problems
     call atomic_jovian(nstat, dmtrx, cgmat, cemat, somat, cumat)

     call atomic_make_amtrx(norbs, cemat, somat, elocs, amtrx)
!     elocs(1:4) = lamb/2.0D0 - 5.2923489070; elocs(5:6) = - lamb - 5.2923489070
!     amtrx = cgmat; elocs(1:4) = lamb/2.0D0; elocs(5:6) = - lamb
!     call zmat_zgemm2(norbs, dmtrx, cgmat, amtrx)
     call atomic_dump_elocs(norbs, elocs)
     call atomic_dump_amtrx(norbs, amtrx)

     print *, "make matrix in imaginary basis ........."
! make angular momentum matrix in imag-axis single particle basis
     call atomic_make_lamat(nband, nspin, norbs, lxmat_i, lymat_i, lzmat_i)
     call atomic_make_samat(nband, nspin, norbs, sxmat_i, symat_i, szmat_i)    
     jxmat_i = lxmat_i + sxmat_i; jymat_i = lymat_i + symat_i; jzmat_i = lzmat_i + szmat_i
! the point group operator C4 
     call atomic_make_c4(nband, nspin, norbs, c4mat_i)

     print *, "end."

     print *, "make matrix in real basis ........."

! make angular momentum matrix in real-axis single particle basis
     call atomic_tran_lamat(norbs, dmtrx, lxmat_i, lymat_i, lzmat_i, lxmat_r, lymat_r, lzmat_r)
     call atomic_tran_samat(norbs, dmtrx, sxmat_i, symat_i, szmat_i, sxmat_r, symat_r, szmat_r)
     call atomic_tran_c4mat(norbs, dmtrx, c4mat_i, c4mat_r) 
     call atomic_tran_jamat(norbs, dmtrx, jxmat_i, jymat_i, jzmat_i, jxmat_r, jymat_r, jzmat_r)
     call atomic_dump_c4mat2(norbs, c4mat_r, 'test_c4_2.dat')
     call atomic_dump_c4mat2(norbs, lxmat_r, 'test_lxr.dat')
     call atomic_dump_c4mat2(norbs, lymat_r, 'test_lyr_2.dat')
     call atomic_dump_c4mat2(norbs, lzmat_r, 'test_lzr_2.dat')

! make c2 sigma_v sigma_d operator in real space
     call atomic_make_c2r(nband, nspin, norbs, c2mat_r)
     call atomic_make_sgmvr(nband, nspin, norbs, sgmvmat_r)
     call atomic_make_sgmdr(nband, nspin, norbs, sgmdmat_r)

     print *, "end."

     print *, "transform matrix in original basis to natural basis ........."

! transformatiom from orginal to natural single particle basis
     call atomic_tran_lamat(norbs, amtrx, lxmat_r, lymat_r, lzmat_r, lxmat_t, lymat_t, lzmat_t)
     call atomic_tran_samat(norbs, amtrx, sxmat_r, symat_r, szmat_r, sxmat_t, symat_t, szmat_t)
     call atomic_tran_c4mat(norbs, amtrx, c4mat_r, c4mat_t) 
     call atomic_tran_c4mat(norbs, amtrx, c2mat_r, c2mat_t) 
     call atomic_tran_c4mat(norbs, amtrx, sgmvmat_r, sgmvmat_t) 
     call atomic_tran_c4mat(norbs, amtrx, sgmdmat_r, sgmdmat_t) 
     call atomic_tran_jamat(norbs, amtrx, jxmat_r, jymat_r, jzmat_r, jxmat_t, jymat_t, jzmat_t)
!     call atomic_dump_c4mat2(norbs, c4mat_t, 'test_c4_nt.dat')

     print *, "end."

     print *, "make Fock basis sheet ........."

! construct Fock basis sheet (just follow nhtong's note)
     call atomic_make_basis( norbs, ncfgs, ntots, nstat, basis, invcd, invsn )

     print *, "end."

     print *, "make matrix in many body space ........."

! build density operator matrix in many body basis
     call atomic_make_npmat(norbs, ncfgs, invcd, npmat)
     
! build single occupied number operator in many body basis
     call atomic_make_sgnpmat(norbs, ncfgs, invcd, sgnpmat)

! build orbital angular momentum operator matrix in many body basis
     call atomic_make_sp2np(norbs, ncfgs, basis, invsn, invcd, lxmat_t, lxmat_n)
     call atomic_make_sp2np(norbs, ncfgs, basis, invsn, invcd, lymat_t, lymat_n)
     call atomic_make_sp2np(norbs, ncfgs, basis, invsn, invcd, lzmat_t, lzmat_n)
     call atomic_dump_c4mat2(ncfgs,   lxmat_n, 'atom.lx.dat')
     call atomic_dump_c4mat2(ncfgs,   lymat_n, 'atom.ly.dat')
     call atomic_dump_c4mat2(ncfgs,   lzmat_n, 'atom.lz.dat')
! build spin angular momentum operator matrix in many body basis
     call atomic_make_sp2np(norbs, ncfgs, basis, invsn, invcd, sxmat_t, sxmat_n)
     call atomic_make_sp2np(norbs, ncfgs, basis, invsn, invcd, symat_t, symat_n)
     call atomic_make_sp2np(norbs, ncfgs, basis, invsn, invcd, szmat_t, szmat_n)
     call atomic_dump_c4mat2(ncfgs,   sxmat_n, 'atom.sx.dat')
     call atomic_dump_c4mat2(ncfgs,   symat_n, 'atom.sy.dat')
     call atomic_dump_c4mat2(ncfgs,   szmat_n, 'atom.sz.dat')

! build c4 operator matrix in many body basis
!     call atomic_make_sp2np(norbs, ncfgs, basis, invsn, invcd, c4mat_t, c4mat_n)
     call gw_make_newUI2(norbs, ncfgs, 5, 7,   c4mat_t, basis, invcd, invsn,   c4mat_n)
     call gw_make_newUI2(norbs, ncfgs, 5, 7,   c2mat_t, basis, invcd, invsn,   c2mat_n)
     call gw_make_newUI2(norbs, ncfgs, 5, 7, sgmvmat_t, basis, invcd, invsn, sgmvmat_n)
     call gw_make_newUI2(norbs, ncfgs, 5, 7, sgmdmat_t, basis, invcd, invsn, sgmdmat_n)
     call atomic_dump_c4mat2(ncfgs,   c4mat_n, 'atom.c4mat.dat')
     call atomic_dump_c4mat2(ncfgs,   c2mat_n, 'atom.c2mat.dat')
     call atomic_dump_c4mat2(ncfgs, sgmvmat_n, 'atom.sgmvmat.dat')
     call atomic_dump_c4mat2(ncfgs, sgmdmat_n, 'atom.sgmdmat.dat')

!     print *, "makeing block-diagonal form C4v ........."
!
!! transform actually block-diagonal matrix into block-diagonal matrix
!     call atomic_tran_blockdiag(ncfgs, c4mat_n, c4mat_bd, c4v_umt, c4degs)
!     call atomic_tran_c4mat(ncfgs, c4v_umt, c2mat_n, c2mat_bd) 
!     call atomic_tran_c4mat(ncfgs, c4v_umt, sgmvmat_n, sgmvmat_bd) 
!     call atomic_tran_c4mat(ncfgs, c4v_umt, sgmdmat_n, sgmdmat_bd) 
!     call atomic_checkirrep(ncfgs, c4mat_bd, c2mat_bd, sgmvmat_bd, sgmdmat_bd, c4degs, irrep_flag, irrep_type, num_rdu)
!     call atomic_dump_c4mat2(ncfgs,   c4mat_bd,    'atom.c4mat_bd.dat')
!     call atomic_dump_c4mat2(ncfgs,   c2mat_bd,   'atom.c2mat_bd.dat')
!     call atomic_dump_c4mat2(ncfgs,   sgmvmat_bd, 'atom.sgmvmat_bd.dat')
!     call atomic_dump_c4mat2(ncfgs,   sgmdmat_bd, 'atom.sgmdmat_bd.dat')
!     call atomic_dump_irrep_flag(ncfgs, c4degs, irrep_flag, irrep_type)

! build total angular momentum operator matrix in many body basis
     call atomic_make_sp2np(norbs, ncfgs, basis, invsn, invcd, jxmat_t, jxmat_n)
     call atomic_make_sp2np(norbs, ncfgs, basis, invsn, invcd, jymat_t, jymat_n)
     call atomic_make_sp2np(norbs, ncfgs, basis, invsn, invcd, jzmat_t, jzmat_n)

! note jjmat and jzmat_n should commute with Hamiltonian matrix
     call atomic_make_l2mat(ncfgs, lxmat_n, lymat_n, lzmat_n, llmat_n)
     call atomic_make_s2mat(ncfgs, sxmat_n, symat_n, szmat_n, ssmat_n)
     call atomic_make_j2mat(ncfgs, jxmat_n, jymat_n, jzmat_n, jjmat_n)
     call atomic_dump_tamat(ncfgs, llmat_n, ssmat_n, jjmat_n, lzmat_n, szmat_n, jzmat_n)

! construct atomic Hamiltonian matrix in {j2, jz} many body basis
     call atomic_tran_cumat(norbs, amtrx, cumat, cumat_t)
     call atomic_make_Hmtrx(norbs, ncfgs, basis, invcd, invsn, elocs, cumat_t, hmat)

     print *, "end."

     print *, "making block-diagonal form C4v ........."

! transform actually block-diagonal matrix into block-diagonal matrix
     call atomic_tran_blockdiag(ncfgs, c4mat_n,  c4mat_bd, c4v_umt, c4v_umt_1d, c4degs, num_2d)

     print *, "end."

     print *, "transferring every operators of C4v into block-diagonal form ........."

     call atomic_dump_c4mat2(ncfgs,    c4v_umt,  'atom.c4v_umt.dat')
     call atomic_tran_c4mat( ncfgs,    c4v_umt,  c2mat_n,   c2mat_bd) 
     call atomic_tran_c4mat( ncfgs,    c4v_umt,  sgmvmat_n, sgmvmat_bd) 
     call atomic_tran_c4mat( ncfgs,    c4v_umt,  sgmdmat_n, sgmdmat_bd) 
     call atomic_tran_c4mat( ncfgs,    c4v_umt,  ssmat_n,   ssmat_bd) 

     print *, "end."

     print *, "checking the irreducibile representation of block-diagonal of C4v ........."
     call atomic_checkirrep( ncfgs,    c4mat_bd, c2mat_bd, sgmvmat_bd, sgmdmat_bd, &
         c4degs, irrep_flag, irrep_type, num_rdu, num_irrep)
     print *, "end."

     print *, "dumping the block-diagonal matrix ........."
     call atomic_dump_c4mat2(ncfgs,    ssmat_bd,  'atom.ssmat_bd.dat')
     call atomic_dump_c4mat2(ncfgs,    c4mat_bd,   'atom.c4mat_bd.dat')
     call atomic_dump_c4mat2(ncfgs,    c2mat_bd,   'atom.c2mat_bd.dat')
     call atomic_dump_c4mat2(ncfgs,    sgmvmat_bd, 'atom.sgmvmat_bd.dat')
     call atomic_dump_c4mat2(ncfgs,    sgmdmat_bd, 'atom.sgmdmat_bd.dat')
     call atomic_dump_irrep_flag(ncfgs, c4degs, irrep_flag, irrep_type,"test_irrep_type.dat")
     print *, "end."

     print *, "showing the details of reducible rep.  ........."
     call atomic_take_reduce(ncfgs, norbs, c4v_umt_1d, irrep_flag, irrep_type, invcd, basis,&
         c4mat_bd, c2mat_bd, sgmvmat_bd, sgmdmat_bd, num_rdu)
     print *, "end."

     print *, "making the unitary matrix for reduction(not general, just for C4v) ........."
     call atomic_reduction_unitary(c4v_rdct_umt)
     print *, "end"

     print *, "making irredcible representation with respect to C4v  ........."
     call atomic_irrep_umtrx(ncfgs, irrep_flag, irrep_type, num_irrep, c4degs, &
         c4v_rdct_umt, c4mat_bd, c2mat_bd, sgmvmat_bd, sgmdmat_bd, c4v_umt, c4virrepdegs)
     call atomic_dump_isgnmat(ncfgs, c4virrepdegs, "atom.c4vdegs.dat")
     call atomic_dump_c4mat2(ncfgs,    c4v_umt,    'atom.c4v_umt_irrep.dat')
     print *, "end."

     print *, "checking ........"
     call atomic_tran_c4mat( ncfgs, c4v_umt,  c4mat_n, c4mat_irrep) 
     call atomic_tran_c4mat( ncfgs, c4v_umt,  c2mat_n, c2mat_irrep) 
     call atomic_tran_c4mat( ncfgs, c4v_umt,  sgmvmat_n, sgmvmat_irrep) 
     call atomic_tran_c4mat( ncfgs, c4v_umt,  sgmdmat_n, sgmdmat_irrep) 
     call atomic_dump_c4mat2(ncfgs, c4mat_irrep,   'atom.c4mat_irrep.dat')
     call atomic_dump_c4mat2(ncfgs, c2mat_irrep,   'atom.c2mat_irrep.dat')
     call atomic_dump_c4mat2(ncfgs, sgmvmat_irrep, 'atom.sgmvmat_irrep.dat')
     call atomic_dump_c4mat2(ncfgs, sgmdmat_irrep, 'atom.sgmdmat_irrep.dat')
     
!     print *, "transform matrix in Fock basis into irrep. of C4v ........."
     call atomic_tran_c4mat(ncfgs,  c4v_umt, szmat_n, szmat_irrep)
     call atomic_tran_c4mat(ncfgs,  c4v_umt, ssmat_n, ssmat_irrep)
     call atomic_tran_c4mat(ncfgs,  c4v_umt, npmat,   npmat_irrep)
     call atomic_tran_c4mat(ncfgs,  c4v_umt, sgnpmat, sgnpmat_irrep)
     call atomic_tran_c4mat(ncfgs,  c4v_umt, hmat,    hmat_irrep)
     call atomic_dump_c4mat2(ncfgs, szmat_irrep, 'atom.szmat_irrep.dat')
     call atomic_dump_c4mat2(ncfgs, ssmat_irrep, 'atom.ssmat_irrep.dat')
     call atomic_dump_c4mat2(ncfgs, sgnpmat_irrep, 'atom.sgnpmat_irrep.dat')
     call atomic_dump_c4mat2(ncfgs, hmat_irrep, 'atom.hmat_irrep.dat')
     call atomic_checkirrep2( ncfgs, c4v_umt, c4mat_irrep, c2mat_irrep, sgmvmat_irrep, &
         sgmdmat_irrep, num_irrep)
     print *, "end."

     print *, "entering diagonalizing Hamiltonian ........."
!     call atomic_diag_Hmat3(norbs, ncfgs, nstat, jzmat_n, jjmat_n, hmat, jzeig, jjeig, Heigs, Heigv)
!     call atomic_diag_Hmat3(norbs, ncfgs, nstat, npmat, sgnpmat, c4mat_n,szmat_n, ssmat_n, hmat, neigs, c4eig, szeig, sseig, Heigs, Heigv)
!     call atomic_diag_Hmat3(norbs, ncfgs, nstat, npmat, sgnpmat, c4mat_n,szmat_n, &
!         ssmat_n, hmat, neigs, szeig, sseig, Heigs, Heigv)
     call atomic_diag_Hmat4(norbs, ncfgs, nstat, c4virrepdegs, c4v_umt, npmat_irrep, &
         sgnpmat_irrep, szmat_irrep, ssmat_irrep, hmat_irrep, neigs, sgeigs, szeig, &
         sseig, Heigs, Heigv, Hdegs)
!     call atomic_dump_taeig(ncfgs, neigs, Heigs, sseig, szeig)
     call atomic_dump_taeig2(ncfgs, irrep_type, neigs, sgeigs, szeig, sseig, Heigs)
     print *, "end."
!     call atomic_make_phase(ncfgs, sxmat_n, symat_n, Heigv)

! dump eigenvalues to file solver.eigs.out
     print *,"entering dumping eigenvalues of atomic hamiltonian ........."
     call atomic_dump_Heigs(ncfgs, heigs)
     print *,"end."


! dump eigenstates to file solver.evec.out
     print *,"entering dumping eigenfunctions of atomic hamiltonian ........."
     call atomic_dump_Heigv(norbs, ncfgs, heigv, invcd)
     print *,"end."

!     call atomic_make_neigs(ncfgs, npmat, Heigv, neigs)
     print *,"entering making vpms ........."
     call atomic_make_vpmmat2(ncfgs, Hdegs, nvpms, vpmmat)
     print *,"end."
     print *,"-------------------------------"

     print *,"the number of reducible of C4v in Fock basis is :"
     print *,num_rdu,'/',num_2d
     print *,"the number of vpms is :"
     print *,nvpms

     print *,"------------------------------"

! deallocate memory and finalize them
     call atomic_final_array()

! print the ending information for atomic problems
     call atomic_print_footer()

     stop
  end program atomic_main

  subroutine atomic_dump_upara(Uc, Uv, Jz, Js, Jp, lamb)
     use constants
     implicit none

     real(dp), intent(in) :: Uc
     real(dp), intent(in) :: Uv
     real(dp), intent(in) :: Jz
     real(dp), intent(in) :: Js
     real(dp), intent(in) :: Jp
     real(dp), intent(in) :: lamb

     open(mytmp, file="atom.Upara.in")
     write(mytmp, "(F17.10)") Uc
     write(mytmp, "(F17.10)") Uv
     write(mytmp, "(F17.10)") Jz
     write(mytmp, "(F17.10)") Js
     write(mytmp, "(F17.10)") Jp
     write(mytmp, "(F17.10)") - lamb
     close(mytmp)

     return
  end subroutine atomic_dump_upara
