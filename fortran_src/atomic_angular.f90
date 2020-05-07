!>>> make lplus operator matrix in imag-axis single particle basis
  subroutine atomic_make_lpmat(nband, nspin, norbs, lpmat)

     implicit none

! number of bands
     integer, intent(in) :: nband

! number of spins
     integer, intent(in) :: nspin

! number of orbits
     integer, intent(in) :: norbs

! jz matrix in single body basis
     complex(8), intent(out) :: lpmat(norbs, norbs)

! local variables
! loop index over orbits
     integer :: iorb, jorb

! orbit angular and spin angular quantum number
! lz ranged from (l:-l),  sz ranged from (s:-s)
     real(8) :: ang_l, ang_s

! orbit and spin momentum value along z axis
     real(8) :: lzval, szval

! initialize lpmat to be zero
     lpmat = dcmplx(0.0D0, 0.0D0)

! orbital and spin angular quantum number
     ang_l = dble(nband-1) / 2.0D0
     ang_s = dble(nspin-1) / 2.0D0

     do jorb=1,norbs
         lzval = ang_l - dble(   (jorb-1)/nspin)
         szval = ang_s - dble(mod(jorb-1,nspin))

         if (abs(lzval - ang_l) .lt. 1.0D-4) cycle
         !$ iorb = jorb - nint(done*(dtwo*ang_s+done))
         iorb = nint((ang_l-(lzval+1.0D0))*(2.0D0*ang_s+1.0D0)) + nint(ang_s-szval) + 1

         lpmat(iorb, jorb) = dsqrt(ang_l*(ang_l+1.0D0) - lzval*(lzval+1.0D0))
     enddo ! over iorb={1,norbs} loop

     return
  end subroutine atomic_make_lpmat

!>>> make lminus operator matrix in imag-axis single particle basis
  subroutine atomic_make_lmmat(nband, nspin, norbs, lmmat)

     implicit none

! external arguments
! number of bands
     integer, intent(in) :: nband

! number of spins
     integer, intent(in) :: nspin

! number of orbits
     integer, intent(in) :: norbs

! lminus matrix in single body basis
     complex(8), intent(out) :: lmmat(norbs, norbs)

! local variables
! loop index over orbits
     integer :: iorb, jorb

! orbit angular and spin angular quantum number
     real(8) :: ang_l, ang_s

! orbit and spin momentum value along z axis
     real(8) :: lzval, szval

! initialize lpmat to be zero
     lmmat = dcmplx(0.0D0, 0.0D0)

! orbital and spin angular quantum number
     ang_l = dble(nband-1) / 2.0D0
     ang_s = dble(nspin-1) / 2.0D0

     do jorb=1,norbs
         lzval = ang_l - dble(   (jorb-1)/nspin)
         szval = ang_s - dble(mod(jorb-1,nspin))

         if (abs(lzval + ang_l) .lt. 1.0D-4) cycle
         iorb = nint((ang_l-(lzval-1.0D0))*(2.0D0*ang_s+1.0D0)) + nint(ang_s-szval) + 1

         lmmat(iorb, jorb) = dsqrt(ang_l*(ang_l+1.0D0) - lzval*(lzval-1.0D0))
     enddo ! over iorb={1,norbs} loop

     return
  end subroutine atomic_make_lmmat

!>>> make lz operator matrix in imag-axis single particle basis
  subroutine atomic_make_lzmat(nband, nspin, norbs, lzmat)

     implicit none

! number of bands
     integer, intent(in) :: nband

! number of spins
     integer, intent(in) :: nspin

! number of orbits
     integer, intent(in) :: norbs

! lz matrix in single body basis
     complex(8), intent(out) :: lzmat(norbs, norbs)

! local variables
! loop index over orbits
     integer :: iorb, jorb

! orbit angular and spin angular quantum number
     real(8) :: ang_l, ang_s

! orbit and spin momentum value along z axis
     real(8) :: lzval, szval

! initialize lpmat to be zero
     lzmat = dcmplx(0.0D0, 0.0D0)

! orbital and spin angular quantum number
     ang_l = dble(nband-1) / 2.0D0
     ang_s = dble(nspin-1) / 2.0D0

     do jorb=1,norbs
         lzval = ang_l - dble(   (jorb-1)/nspin)
         szval = ang_s - dble(mod(jorb-1,nspin))
         iorb = jorb; lzmat(jorb, jorb) = lzval
     enddo ! over iorb={1,norbs} loop

     return
  end subroutine atomic_make_lzmat

!>>> make splus operator matrix in axis single particle basis
  subroutine atomic_make_spmat(nband, nspin, norbs, spmat)

     implicit none

! number of bands
     integer, intent(in) :: nband

! number of spins
     integer, intent(in) :: nspin

! number of orbits
     integer, intent(in) :: norbs

! splus matrix in single body basis
     complex(8), intent(out) :: spmat(norbs, norbs)

! local variables
! loop index over orbits
     integer :: iorb, jorb

! orbit angular and spin angular quantum number
     real(8) :: ang_l, ang_s

! orbit and spin momentum value along z axis
     real(8) :: lzval, szval

! initialize lpmat to be zero
     spmat = dcmplx(0.0D0, 0.0D0)

! orbital and spin angular quantum number
     ang_l = dble(nband-1) / 2.0D0
     ang_s = dble(nspin-1) / 2.0D0

     do jorb=1,norbs
         lzval = ang_l - dble(   (jorb-1)/nspin)
         szval = ang_s - dble(mod(jorb-1,nspin))

         if (abs(szval - ang_s) .lt. 1.0D-4 ) cycle
         iorb = nint((ang_l-lzval)*(2.0D0*ang_s+1.0D0)) + nint(ang_s-szval-1.0D0) + 1

         spmat(iorb, jorb) = dsqrt(ang_s*(ang_s+1.0D0) - szval*(szval+1.0D0))
     enddo ! over iorb={1,norbs} loop

     return
  end subroutine atomic_make_spmat

!>>> make sminus operator matrix in imag-axis single particle basis
  subroutine atomic_make_smmat(nband, nspin, norbs, smmat)

     implicit none

! number of bands
     integer, intent(in) :: nband

! number of spins
     integer, intent(in) :: nspin

! number of orbits
     integer, intent(in) :: norbs

! sminus matrix in single body basis
     complex(8), intent(out) :: smmat(norbs, norbs)

! local variables
! loop index over orbits
     integer :: iorb, jorb

! orbit angular and spin angular quantum number
     real(8) :: ang_l, ang_s

! orbit and spin momentum value along z axis
     real(8) :: lzval, szval

! initialize smmat to be zero
     smmat = dcmplx(0.0D0, 0.0D0)

! orbital and spin angular quantum number
     ang_l = dble(nband-1) / 2.0D0
     ang_s = dble(nspin-1) / 2.0D0

     do jorb=1,norbs
         lzval = ang_l - dble(   (jorb-1)/nspin)
         szval = ang_s - dble(mod(jorb-1,nspin))

         if (abs(szval + ang_s) .lt. 1.0D-4) cycle
         iorb = nint((ang_l-lzval)*(2.0D0*ang_s+1.0D0)) + nint(ang_s-szval+1.0D0) + 1

         smmat(iorb, jorb) = dsqrt(ang_s*(ang_s+1.0D0) - szval*(szval-1.0D0))
     enddo ! over iorb={1,norbs} loop

     return
  end subroutine atomic_make_smmat

!>>> make sz operator matrix in imag-axis single particle basis
  subroutine atomic_make_szmat(nband, nspin, norbs, szmat)

     implicit none

! number of bands
     integer, intent(in) :: nband

! number of spins
     integer, intent(in) :: nspin

! number of orbits
     integer, intent(in) :: norbs

! jz matrix in single body basis
     complex(8), intent(out) :: szmat(norbs, norbs)

! local variables
! loop index over orbits
     integer :: iorb, jorb

! orbit angular and spin angular quantum number
     real(8) :: ang_l, ang_s

! orbit and spin momentum value along z axis
     real(8) :: lzval, szval

! initialize lpmat to be zero
     szmat = dcmplx(0.0D0, 0.0D0)

! orbital and spin angular quantum number
     ang_l = dble(nband-1) / 2.0D0
     ang_s = dble(nspin-1) / 2.0D0

     do jorb=1,norbs
         lzval = ang_l - dble(   (jorb-1)/nspin)
         szval = ang_s - dble(mod(jorb-1,nspin))
         iorb = jorb; szmat(jorb, jorb) = szval
     enddo ! over iorb={1,norbs} loop

     return
  end subroutine atomic_make_szmat

  subroutine atomic_make_l2mat(ncfgs, lxmat, lymat, lzmat, l2mat)
     implicit none

! number of configurations
     integer, intent(in) :: ncfgs

! lx operator in many particle basis
     complex(8), intent(in) :: lxmat(ncfgs, ncfgs)

! ly operator in many particle basis
     complex(8), intent(in) :: lymat(ncfgs, ncfgs)

! lz operator in many particle basis
     complex(8), intent(in) :: lzmat(ncfgs, ncfgs)

! l2 operator in many particle basis
     complex(8), intent(out) :: l2mat(ncfgs, ncfgs)

     l2mat = dcmplx(0.0D0, 0.0D0)
     l2mat = l2mat + matmul(lxmat, lxmat)
     l2mat = l2mat + matmul(lymat, lymat)
     l2mat = l2mat + matmul(lzmat, lzmat)

     return
  end subroutine atomic_make_l2mat

  subroutine atomic_make_s2mat(ncfgs, sxmat, symat, szmat, s2mat)
     implicit none

! number of configurations
     integer, intent(in) :: ncfgs

! sx operator in many particle basis
     complex(8), intent(in) :: sxmat(ncfgs, ncfgs)

! sy operator in many particle basis
     complex(8), intent(in) :: symat(ncfgs, ncfgs)

! sz operator in many particle basis
     complex(8), intent(in) :: szmat(ncfgs, ncfgs)

! s2 operator in many particle basis
     complex(8), intent(out) :: s2mat(ncfgs, ncfgs)

     s2mat = dcmplx(0.0D0, 0.0D0)
     s2mat = s2mat + matmul(sxmat, sxmat)
     s2mat = s2mat + matmul(symat, symat)
     s2mat = s2mat + matmul(szmat, szmat)

     return
  end subroutine atomic_make_s2mat

  subroutine atomic_make_j2mat(ncfgs, jxmat, jymat, jzmat, j2mat)
     implicit none

! number of configurations
     integer, intent(in) :: ncfgs

! jx operator in many particle basis
     complex(8), intent(in) :: jxmat(ncfgs, ncfgs)

! jy operator in many particle basis
     complex(8), intent(in) :: jymat(ncfgs, ncfgs)

! jz operator in many particle basis
     complex(8), intent(in) :: jzmat(ncfgs, ncfgs)

! j2 operator in many particle basis
     complex(8), intent(out) :: j2mat(ncfgs, ncfgs)

     j2mat = dcmplx(0.0D0, 0.0D0)
     j2mat = j2mat + matmul(jxmat, jxmat)
     j2mat = j2mat + matmul(jymat, jymat)
     j2mat = j2mat + matmul(jzmat, jzmat)

     return
  end subroutine atomic_make_j2mat

  subroutine atomic_make_dmatp(nband, nspin, norbs, dmat)

     implicit none

! external arguments
! number of band
     integer, intent(in) :: nband

! number of spin
     integer, intent(in) :: nspin

! number of orbits
     integer, intent(in) :: norbs

! transfromation from real-axis to imag-axis single particle basis
     complex(8), intent(out) :: dmat(norbs, norbs)

! loop index over nband
     integer :: ibnd
     integer :: jbnd

! auxiliary integer arrays
     integer, allocatable :: ipiv(:)

! allocate memory for ipiv and initialize it
     allocate(ipiv(norbs))
     do ibnd=1,nband
         jbnd = ibnd + nband
         ipiv(ibnd) = 2 * ibnd - 1
         ipiv(jbnd) = 2 * ibnd - 0
     enddo ! over ibnd={1,nband} loop

! initialize the transformation matrix
     dmat = dcmplx(0.0D0, 0.0D0)

! left top block diagonal parts
     dmat(ipiv(1), ipiv(1)) = - dcmplx(1.0D0, 0.0D0)/dsqrt(2.0D0)
     dmat(ipiv(3), ipiv(1)) = + dcmplx(1.0D0, 0.0D0)/dsqrt(2.0D0)
     dmat(ipiv(1), ipiv(2)) = + dcmplx(0.0D0, 1.0D0)/dsqrt(2.0D0)
     dmat(ipiv(3), ipiv(2)) = + dcmplx(0.0D0, 1.0D0)/dsqrt(2.0D0)
     dmat(ipiv(2), ipiv(3)) = + dcmplx(1.0D0, 0.0D0)/dsqrt(1.0D0)

! right bottom block diagonal parts
     dmat(ipiv(1+nband), ipiv(1+nband)) = dmat(ipiv(1), ipiv(1))
     dmat(ipiv(3+nband), ipiv(1+nband)) = dmat(ipiv(3), ipiv(1))
     dmat(ipiv(1+nband), ipiv(2+nband)) = dmat(ipiv(1), ipiv(2))
     dmat(ipiv(3+nband), ipiv(2+nband)) = dmat(ipiv(3), ipiv(2))
     dmat(ipiv(2+nband), ipiv(3+nband)) = dmat(ipiv(2), ipiv(3))

     return
  end subroutine atomic_make_dmatp

  subroutine atomic_make_dmatd(nband, nspin, norbs, dmat)

     implicit none

! external arguments
! number of band
     integer, intent(in) :: nband

! number of spin
     integer, intent(in) :: nspin

! number of orbits
     integer, intent(in) :: norbs

! transfromation from imag-axis single particle basis to real basis
! basis order : 2,1,0,-1,-2 ----> dyz, dxz, dxy, dz2, dx2-y2
! spin order  : ududududud
     complex(8), intent(out) :: dmat(norbs, norbs)

! loop index over nband
     integer :: ibnd
     integer :: jbnd

! auxiliary integer arrays
     integer, allocatable :: ipiv(:)

! allocate memory for ipiv and initialize it
     allocate(ipiv(norbs)); ipiv = 0
     do ibnd=1,nband
         jbnd = ibnd + nband
         ipiv(ibnd) = 2 * ibnd - 1
         ipiv(jbnd) = 2 * ibnd - 0
     enddo ! over ibnd={1,nband} loop

! initialize the transformation matrix
     dmat = dcmplx(0.0D0, 0.0D0)

! left top block diagonal parts
     dmat(ipiv(2), ipiv(1)) = + dcmplx(0.0D0, 1.0D0)/dsqrt(2.0D0)
     dmat(ipiv(4), ipiv(1)) = + dcmplx(0.0D0, 1.0D0)/dsqrt(2.0D0)
     dmat(ipiv(2), ipiv(2)) = - dcmplx(1.0D0, 0.0D0)/dsqrt(2.0D0)
     dmat(ipiv(4), ipiv(2)) = + dcmplx(1.0D0, 0.0D0)/dsqrt(2.0D0)
     dmat(ipiv(1), ipiv(3)) = - dcmplx(0.0D0, 1.0D0)/dsqrt(2.0D0)
     dmat(ipiv(5), ipiv(3)) = + dcmplx(0.0D0, 1.0D0)/dsqrt(2.0D0)
     dmat(ipiv(3), ipiv(4)) = + dcmplx(1.0D0, 0.0D0)/dsqrt(1.0D0)
     dmat(ipiv(1), ipiv(5)) = + dcmplx(1.0D0, 0.0D0)/dsqrt(2.0D0)
     dmat(ipiv(5), ipiv(5)) = + dcmplx(1.0D0, 0.0D0)/dsqrt(2.0D0)

! right bottom block diagonal parts
     dmat(ipiv(2+nband), ipiv(1+nband)) = dmat(ipiv(2), ipiv(1))
     dmat(ipiv(4+nband), ipiv(1+nband)) = dmat(ipiv(4), ipiv(1))
     dmat(ipiv(2+nband), ipiv(2+nband)) = dmat(ipiv(2), ipiv(2))
     dmat(ipiv(4+nband), ipiv(2+nband)) = dmat(ipiv(4), ipiv(2))
     dmat(ipiv(1+nband), ipiv(3+nband)) = dmat(ipiv(1), ipiv(3))
     dmat(ipiv(5+nband), ipiv(3+nband)) = dmat(ipiv(5), ipiv(3))
     dmat(ipiv(3+nband), ipiv(4+nband)) = dmat(ipiv(3), ipiv(4))
     dmat(ipiv(1+nband), ipiv(5+nband)) = dmat(ipiv(1), ipiv(5))
     dmat(ipiv(5+nband), ipiv(5+nband)) = dmat(ipiv(5), ipiv(5))

! deallocate auxiliary arrays
     if (allocated(ipiv)) deallocate(ipiv)

     return
  end subroutine atomic_make_dmatd

! make Orbital Angular Momentum matrix in imag-axis single particle basis
  subroutine atomic_make_lamat(nband, nspin, norbs, lxmat, lymat, lzmat)

     implicit none

! number of bands
     integer, intent(in) :: nband

! number of spins
     integer, intent(in) :: nspin

! number of orbits
     integer, intent(in) :: norbs

! lx matrix in real-axis single body basis
     complex(8), intent(out) :: lxmat(norbs, norbs)

! ly matrix in real-axis single body basis
     complex(8), intent(out) :: lymat(norbs, norbs)

! lz matrix in real-axis single body basis
     complex(8), intent(out) :: lzmat(norbs, norbs)

! local variables
! status while allocating memory
     integer :: istat

! auxiliary complex(8) allocatable arrays
     complex(8), allocatable :: lpmat(:, :)
     complex(8), allocatable :: lmmat(:, :)

! allocate memory for allocatable arrays
     allocate(lpmat(norbs, norbs), lmmat(norbs, norbs), stat=istat)
     if (istat .ne. 0) stop "allocate memory error in atomic_make_LAmat"

! initialize the allocated arrays
     lpmat = dcmplx(0.0D0, 0.0D0)
     lmmat = dcmplx(0.0D0, 0.0D0)

! construct Orbital Angular momentum operator along xy axis
     call atomic_make_lpmat(nband, nspin, norbs, lpmat)
     call atomic_make_lmmat(nband, nspin, norbs, lmmat)

     lxmat = (lpmat + lmmat) / dcmplx(2.0D0, 0.0D0)
     lymat = (lpmat - lmmat) / dcmplx(0.0D0, 2.0D0)

! construct Orbital Angular momentum operator along z  axis
     call atomic_make_lzmat(nband, nspin, norbs, lzmat)

! deallocate memory for allocated auxiliary arrays
     if (allocated(lpmat)) deallocate(lpmat)
     if (allocated(lmmat)) deallocate(lmmat)

     return
  end subroutine atomic_make_lamat

!>>> transfrom Spin Angular matrix from imag-axis (orginal) to real-axis (natural) basis
  subroutine atomic_tran_samat(norbs, dmat, sxmat, symat, szmat, sxmat_t, symat_t, szmat_t)

     implicit none

! external arguments
! number of orbits
     integer, intent(in) :: norbs

! transformation matrix from imag-axis(orginal) to real-axis(natural) basis
     complex(8), intent(in) :: dmat(norbs, norbs)

! sx matrix in image-axis (orginal) single particle basis
     complex(8), intent(in) :: sxmat(norbs, norbs)

! sy matrix in imag-axis (orginal) single particle basis
     complex(8), intent(in) :: symat(norbs, norbs)

! sz matrix in imag-axis (orginal) single particle basis
     complex(8), intent(in) :: szmat(norbs, norbs)

! sx matrix in real-axis (natural) single particle basis
     complex(8), intent(out) :: sxmat_t(norbs, norbs)

! sy matrix in real-axis (natural) single particle basis
     complex(8), intent(out) :: symat_t(norbs, norbs)

! sz matrix in real-axis (natural) single particle basis
     complex(8), intent(out) :: szmat_t(norbs, norbs)

! unitary transformation for sx operator
     call zmat_unit2(norbs, dmat, sxmat, sxmat_t)

! unitary transformation for sy operator
     call zmat_unit2(norbs, dmat, symat, symat_t)

! unitary transformation for sz operator
     call zmat_unit2(norbs, dmat, szmat, szmat_t)

     return
  end subroutine atomic_tran_samat

!>>> transfrom c4 matrix from imag-axis (orginal) to real-axis (natural) basis
  subroutine atomic_tran_c4mat(norbs, dmat, c4mat, c4mat_t)

     implicit none

! external arguments
! number of orbits
     integer, intent(in) :: norbs

! transformation matrix from imag-axis(orginal) to real-axis(natural) basis
     complex(8), intent(in) :: dmat(norbs, norbs)

! c4 matrix in imag-axis (orginal) single particle basis
     complex(8), intent(in) :: c4mat(norbs, norbs)

! c4 matrix in real-axis (natural) single particle basis
     complex(8), intent(out) :: c4mat_t(norbs, norbs)

! unitary transformation for c4 operator
     call zmat_unit2(norbs, dmat, c4mat, c4mat_t)

     return
  end subroutine atomic_tran_c4mat

!>>> transfrom Orbital Angular matrix from imag-axis (orginal) to real-axis (natural) basis
  subroutine atomic_tran_lamat(norbs, dmat, lxmat, lymat, lzmat, lxmat_t, lymat_t, lzmat_t)

     implicit none

! external arguments
! number of orbits
     integer, intent(in) :: norbs

! transformation matrix from imag-axis(orginal) to real-axis(natural) basis
     complex(8), intent(in) :: dmat(norbs, norbs)

! lx matrix in image-axis (orginal) single particle basis
     complex(8), intent(in) :: lxmat(norbs, norbs)

! ly matrix in imag-axis (orginal) single particle basis
     complex(8), intent(in) :: lymat(norbs, norbs)

! lz matrix in imag-axis (orginal) single particle basis
     complex(8), intent(in) :: lzmat(norbs, norbs)

! lx matrix in real-axis (natural) single particle basis
     complex(8), intent(out) :: lxmat_t(norbs, norbs)

! ly matrix in real-axis (natural) single particle basis
     complex(8), intent(out) :: lymat_t(norbs, norbs)

! lz matrix in real-axis (natural) single particle basis
     complex(8), intent(out) :: lzmat_t(norbs, norbs)

! unitary transformation for sx operator
     call zmat_unit2(norbs, dmat, lxmat, lxmat_t)

! unitary transformation for sy operator
     call zmat_unit2(norbs, dmat, lymat, lymat_t)

! unitary transformation for sz operator
     call zmat_unit2(norbs, dmat, lzmat, lzmat_t)

     return
  end subroutine atomic_tran_lamat

  subroutine atomic_tran_jamat(norbs, dmat, jxmat, jymat, jzmat, jxmat_t, jymat_t, jzmat_t)

     implicit none

! external arguments
! number of orbits
     integer, intent(in) :: norbs

! transformation matrix from imag-axis(orginal) to real-axis(natural) basis
     complex(8), intent(in) :: dmat(norbs, norbs)

! lx matrix in image-axis (orginal) single particle basis
     complex(8), intent(in) :: jxmat(norbs, norbs)

! ly matrix in imag-axis (orginal) single particle basis
     complex(8), intent(in) :: jymat(norbs, norbs)

! lz matrix in imag-axis (orginal) single particle basis
     complex(8), intent(in) :: jzmat(norbs, norbs)

! lx matrix in real-axis (natural) single particle basis
     complex(8), intent(out) :: jxmat_t(norbs, norbs)

! ly matrix in real-axis (natural) single particle basis
     complex(8), intent(out) :: jymat_t(norbs, norbs)

! lz matrix in real-axis (natural) single particle basis
     complex(8), intent(out) :: jzmat_t(norbs, norbs)

! unitary transformation for sx operator
     call zmat_unit2(norbs, dmat, jxmat, jxmat_t)

! unitary transformation for sy operator
     call zmat_unit2(norbs, dmat, jymat, jymat_t)

! unitary transformation for sz operator
     call zmat_unit2(norbs, dmat, jzmat, jzmat_t)

     return
  end subroutine atomic_tran_jamat

!>>> make C4 operator representation in imag-axis single particle basis
  subroutine atomic_make_c4(nband, nspin, norbs, c4mat)

     implicit none

! external arguments
! number of bands
     integer, intent(in) :: nband

! number of spins
     integer, intent(in) :: nspin

! number of orbits
     integer, intent(in) :: norbs

! c4 matrix in real-axis single body basis
     complex(8), intent(out) :: c4mat(norbs, norbs)

! local variables
! status while allocating memory
     integer :: istat
! local integer
     integer :: i,j,k,l,ibase,jbase

! auxiliary complex(8) variables
     complex(8), allocatable :: spmat(:, :)
     complex(8), allocatable :: smmat(:, :)

! allocate memory for allocatable arrays
     allocate(spmat(nband, nband), smmat(nspin, nspin), stat=istat)
     if (istat .ne. 0) stop "allocate memory error in atomic_make_c4mat"

! initialize the allocated arrays
     spmat = dcmplx(0.0D0, 0.0D0)
     smmat = dcmplx(0.0D0, 0.0D0)
     c4mat = dcmplx(0.0D0, 0.0D0)

! check whether d orbitals or not
     if (nband .ne. 5 .or. nspin .ne. 2) stop "other type is not support yet!"

! construct c4 matrix
! c4 in space 
     spmat(1,1)=dcmplx(-1.0D0,  0.0D0)
     spmat(2,2)=dcmplx( 0.0D0, -1.0D0)
     spmat(3,3)=dcmplx( 1.0D0,  0.0D0)
     spmat(4,4)=dcmplx( 0.0D0,  1.0D0)
     spmat(5,5)=dcmplx(-1.0D0,  0.0D0)
! c4 in spin (it is unit when no soc)
     smmat(1,1)=dcmplx( 1.0D0,  0.0D0)
     smmat(2,2)=dcmplx( 1.0D0,  0.0D0)
! c4 in joint space of spin and orbital
     do i = 1, nband
         ibase = (i-1)*nspin
         do j = 1,nband
             jbase = (j-1)*nspin
             do k = 1, nspin
                 do l = 1, nspin
                     c4mat(ibase+k,jbase+l) = spmat(i,j)*smmat(k,l) 
                 enddo
             enddo
         enddo
     enddo

! deallocate memory for allocated arrays
     if (allocated(spmat)) deallocate(spmat)
     if (allocated(smmat)) deallocate(smmat)

     return
  end subroutine atomic_make_c4

!>>> make C2 operator representation in real-axis single particle basis
  subroutine atomic_make_c2r(nband, nspin, norbs, c2mat)

     implicit none

! external arguments
! number of bands
     integer, intent(in) :: nband

! number of spins
     integer, intent(in) :: nspin

! number of orbits
     integer, intent(in) :: norbs

! c4 matrix in real-axis single body basis
     complex(8), intent(out) :: c2mat(norbs, norbs)

! local variables
! status while allocating memory
     integer :: istat
! local integer
     integer :: i,j,k,l,ibase,jbase

! auxiliary complex(8) variables
     complex(8), allocatable :: spmat(:, :)
     complex(8), allocatable :: smmat(:, :)

! allocate memory for allocatable arrays
     allocate(spmat(nband, nband), smmat(nspin, nspin), stat=istat)
     if (istat .ne. 0) stop "allocate memory error in atomic_make_c2mat"

! initialize the allocated arrays
     spmat = dcmplx(0.0D0, 0.0D0)
     smmat = dcmplx(0.0D0, 0.0D0)
     c2mat = dcmplx(0.0D0, 0.0D0)

! check whether d orbitals or not
     if (nband .ne. 5 .or. nspin .ne. 2) stop "other type is not support yet!"

! construct c2 matrix
! c2 in space 
     spmat(1,1)=dcmplx(-1.0D0,  0.0D0)
     spmat(2,2)=dcmplx(-1.0D0,  0.0D0)
     spmat(3,3)=dcmplx( 1.0D0,  0.0D0)
     spmat(4,4)=dcmplx( 1.0D0,  0.0D0)
     spmat(5,5)=dcmplx( 1.0D0,  0.0D0)
! c4 in spin (it is unit when no soc)
     smmat(1,1)=dcmplx( 1.0D0,  0.0D0)
     smmat(2,2)=dcmplx( 1.0D0,  0.0D0)
! c4 in joint space of spin and orbital
     do i = 1, nband
         ibase = (i-1)*nspin
         do j = 1,nband
             jbase = (j-1)*nspin
             do k = 1, nspin
                 do l = 1, nspin
                     c2mat(ibase+k,jbase+l) = spmat(i,j)*smmat(k,l) 
                 enddo
             enddo
         enddo
     enddo

! deallocate memory for allocated arrays
     if (allocated(spmat)) deallocate(spmat)
     if (allocated(smmat)) deallocate(smmat)

     return
  end subroutine atomic_make_c2r

!>>> make sigma_d operator representation in real-axis single particle basis
  subroutine atomic_make_sgmdr(nband, nspin, norbs, sgmdmat)

     implicit none

! external arguments
! number of bands
     integer, intent(in) :: nband

! number of spins
     integer, intent(in) :: nspin

! number of orbits
     integer, intent(in) :: norbs

! c4 matrix in real-axis single body basis
     complex(8), intent(out) :: sgmdmat(norbs, norbs)

! local variables
! status while allocating memory
     integer :: istat
! local integer
     integer :: i,j,k,l,ibase,jbase

! auxiliary complex(8) variables
     complex(8), allocatable :: spmat(:, :)
     complex(8), allocatable :: smmat(:, :)

! allocate memory for allocatable arrays
     allocate(spmat(nband, nband), smmat(nspin, nspin), stat=istat)
     if (istat .ne. 0) stop "allocate memory error in atomic_make_c2mat"

! initialize the allocated arrays
     spmat   = dcmplx(0.0D0, 0.0D0)
     smmat   = dcmplx(0.0D0, 0.0D0)
     sgmdmat = dcmplx(0.0D0, 0.0D0)

! check whether d orbitals or not
     if (nband .ne. 5 .or. nspin .ne. 2) stop "other type is not support yet!"

! construct c2 matrix
! c2 in space 
     spmat(1,2)=dcmplx( 1.0D0,  0.0D0)
     spmat(2,1)=dcmplx( 1.0D0,  0.0D0)
     spmat(3,3)=dcmplx( 1.0D0,  0.0D0)
     spmat(4,4)=dcmplx( 1.0D0,  0.0D0)
     spmat(5,5)=dcmplx(-1.0D0,  0.0D0)
! c4 in spin (it is unit when no soc)
     smmat(1,1)=dcmplx( 1.0D0,  0.0D0)
     smmat(2,2)=dcmplx( 1.0D0,  0.0D0)
! c4 in joint space of spin and orbital
     do i = 1, nband
         ibase = (i-1)*nspin
         do j = 1,nband
             jbase = (j-1)*nspin
             do k = 1, nspin
                 do l = 1, nspin
                     sgmdmat(ibase+k,jbase+l) = spmat(i,j)*smmat(k,l) 
                 enddo
             enddo
         enddo
     enddo

! deallocate memory for allocated arrays
     if (allocated(spmat)) deallocate(spmat)
     if (allocated(smmat)) deallocate(smmat)

     return
  end subroutine atomic_make_sgmdr

!>>> make sigma_d operator representation in real-axis single particle basis
  subroutine atomic_make_sgmvr(nband, nspin, norbs, sgmvmat)

     implicit none

! external arguments
! number of bands
     integer, intent(in) :: nband

! number of spins
     integer, intent(in) :: nspin

! number of orbits
     integer, intent(in) :: norbs

! c4 matrix in real-axis single body basis
     complex(8), intent(out) :: sgmvmat(norbs, norbs)

! local variables
! status while allocating memory
     integer :: istat
! local integer
     integer :: i,j,k,l,ibase,jbase

! auxiliary complex(8) variables
     complex(8), allocatable :: spmat(:, :)
     complex(8), allocatable :: smmat(:, :)

! allocate memory for allocatable arrays
     allocate(spmat(nband, nband), smmat(nspin, nspin), stat=istat)
     if (istat .ne. 0) stop "allocate memory error in atomic_make_c2mat"

! initialize the allocated arrays
     spmat   = dcmplx(0.0D0, 0.0D0)
     smmat   = dcmplx(0.0D0, 0.0D0)
     sgmvmat = dcmplx(0.0D0, 0.0D0)

! check whether d orbitals or not
     if (nband .ne. 5 .or. nspin .ne. 2) stop "other type is not support yet!"

! construct c2 matrix
! c2 in space 
     spmat(1,1)=dcmplx( 1.0D0,  0.0D0)
     spmat(2,2)=dcmplx(-1.0D0,  0.0D0)
     spmat(3,3)=dcmplx(-1.0D0,  0.0D0)
     spmat(4,4)=dcmplx( 1.0D0,  0.0D0)
     spmat(5,5)=dcmplx( 1.0D0,  0.0D0)
! c4 in spin (it is unit when no soc)
     smmat(1,1)=dcmplx( 1.0D0,  0.0D0)
     smmat(2,2)=dcmplx( 1.0D0,  0.0D0)
! c4 in joint space of spin and orbital
     do i = 1, nband
         ibase = (i-1)*nspin
         do j = 1,nband
             jbase = (j-1)*nspin
             do k = 1, nspin
                 do l = 1, nspin
                     sgmvmat(ibase+k,jbase+l) = spmat(i,j)*smmat(k,l) 
                 enddo
             enddo
         enddo
     enddo

! deallocate memory for allocated arrays
     if (allocated(spmat)) deallocate(spmat)
     if (allocated(smmat)) deallocate(smmat)

     return
  end subroutine atomic_make_sgmvr

!>>> make Spin Angular matrix (sx, sy, sz) in imag-axis single particle basis
  subroutine atomic_make_samat(nband, nspin, norbs, sxmat, symat, szmat)

     implicit none

! external arguments
! number of bands
     integer, intent(in) :: nband

! number of spins
     integer, intent(in) :: nspin

! number of orbits
     integer, intent(in) :: norbs

! sx matrix in real-axis single body basis
     complex(8), intent(out) :: sxmat(norbs, norbs)

! sy matrix in real-axis single body basis
     complex(8), intent(out) :: symat(norbs, norbs)

! sz matrix in real-axis single body basis
     complex(8), intent(out) :: szmat(norbs, norbs)

! local variables
! status while allocating memory
     integer :: istat

! auxiliary complex(8) variables
     complex(8), allocatable :: spmat(:, :)
     complex(8), allocatable :: smmat(:, :)

! allocate memory for allocatable arrays
     allocate(spmat(norbs, norbs), smmat(norbs, norbs), stat=istat)
     if (istat .ne. 0) stop "allocate memory error in atomic_make_SAmat"

! initialize the allocated arrays
     spmat = dcmplx(0.0D0, 0.0D0)
     smmat = dcmplx(0.0D0, 0.0D0)

! construct Spin Angular momentum operator
     call atomic_make_spmat(nband, nspin, norbs, spmat)
     call atomic_make_smmat(nband, nspin, norbs, smmat)
     sxmat = (spmat + smmat) / dcmplx(2.0D0, 0.0D0)
     symat = (spmat - smmat) / dcmplx(0.0D0, 2.0D0)

! construct Spin Angular momentum operator along z-axis
     call atomic_make_szmat(nband, nspin, norbs, szmat)

! deallocate memory for allocated arrays
     if (allocated(spmat)) deallocate(spmat)
     if (allocated(smmat)) deallocate(smmat)

     return
  end subroutine atomic_make_samat

!--------------------three band for example------------------------!
!   norbs    ang_l    ang_s    lzval    szval    ang_j    jzval    !
!     1        1       1/2       1       1/2      3/2      3/2     !
!     2        1       1/2       1      -1/2      3/2      1/2     !
!     3        1       1/2       0       1/2      3/2     -1/2     !
!     4        1       1/2       0      -1/2      3/2     -3/2     !
!     5        1       1/2      -1       1/2      1/2      1/2     !
!     6        1       1/2      -1      -1/2      1/2     -1/2     !
!------------------------------------------------------------------!
!     1        1     1.0000000000                                  !
!     2        2     0.5773502692                                  !
!     2        5     0.8164965809                                  !
!     3        2     0.8164965809                                  !
!     3        5    -0.5773502692                                  !
!     4        3     0.8164965809                                  !
!     4        6     0.5773502692                                  !
!     5        3     0.5773502692                                  !
!     5        6    -0.8164965809                                  !
!     6        4     1.0000000000                                  !
!------------------------------------------------------------------!
  subroutine atomic_make_cgmat(nband, nspin, norbs, cgmat)
     use constants

     implicit none

! external arguments
! clebsch-gordan coefficient
     real(8), external :: clebsch_gordan

! number of band
     integer, intent(in) :: nband

! number of spin
     integer, intent(in) :: nspin

! number of orbits
     integer, intent(in) :: norbs

! transformation matrix from {lz, sz} to {j2, jz} basis
     complex(dp), intent(out) :: cgmat(norbs, norbs)

! local variables
! loop index over orbits
     integer :: iorb, jorb

! degeneracy of the two subset in {j2, jz} basis
     integer :: ndeg, mdeg

! orbit angular, spin and total angular momentum quantum number
     real(dp) :: ang_l, ang_s, ang_j

! orbit and spin and total angular momentum value along z axis
     real(dp) :: lzval, szval, jzval

! initialize cgmat to be zero
     cgmat = czero

! orbital and spin angular quantum number
     ang_l = dble(nband-1) / dtwo
     ang_s = dble(nspin-1) / dtwo

! degeneracy of the two total angular momentum quantum number
     ndeg  = nint(dtwo*(ang_l + ang_s) + done)
     mdeg  = nint(dtwo*(ang_l - ang_s) + done)
     if ((ndeg + mdeg) /= norbs) stop "serve error happened in atomic_make_cgmat"
     
     do jorb=1,norbs
         if (jorb .le. ndeg) then
             ang_j = ang_l + ang_s
             jzval = ang_j - dble(jorb - 1)
         else
             ang_j = ang_l - ang_s
             jzval = ang_j - dble(jorb - 1 - ndeg)
         endif ! back if (jorb .le. ndeg0 ) block

         do iorb=1,norbs
             lzval = ang_l - dble(   (iorb-1)/nspin)
             szval = ang_s - dble(mod(iorb-1,nspin))
             if (abs(lzval + szval - jzval) .gt. eps3) cycle
             cgmat(iorb, jorb) = clebsch_gordan(ang_l, ang_s, ang_j, lzval, szval, jzval)
         enddo ! over iorb={1,norbs} loop
     enddo ! over jorb={1,norbs} loop

     return
  end subroutine atomic_make_cgmat

! Operator transformation from single particle basis to many body basis
  subroutine atomic_make_sp2np(norbs, ncfgs, basis, invsn, invcd, zmat_s, zmat_n)
     implicit none

! number of orbits
     integer, intent(in) :: norbs

! number of configurations
     integer, intent(in) :: ncfgs

! decimal represented basis
     integer, intent(in) :: basis(ncfgs)

! serial number of a decimal represented basis
     integer, intent(in) :: invsn(0:2**norbs-1)

! binary represented basis
     integer, intent(in) :: invcd(norbs, ncfgs)

! single particle operator in single particle basis
     complex(8), intent(in) :: zmat_s(norbs, norbs)

! single particle operator in many particle basis
     complex(8), intent(out) :: zmat_n(ncfgs, ncfgs)

! local variables
! loop index over orbits
     integer :: iorb
     integer :: alpha
     integer :: betta

! loop index over configurations
     integer :: ibas
     integer :: jbas

! sign change due to anti-commute relation
     integer :: isgn

! auxiliary integer variables
     integer :: jold, jnew
     integer :: code(norbs)

! auxiliary real(dp) variables
     real(8) :: dsgn

! initialize zmat_n to be zero
     zmat_n = dcmplx(0.0d0, 0.0d0)

! main loop over many body basis
     jbasloop: do jbas=1,ncfgs
         alphaloop: do alpha=1,norbs
         bettaloop: do betta=1,norbs

             ! initialize some variables
             if (abs(zmat_s(alpha, betta)) .lt. 1.0D-6) cycle
             isgn = 0; jold = basis(jbas); code = invcd(:, jbas)

             ! simulate one eliminate operator
             if (code(betta) .eq. 1) then
                 do iorb=1,betta-1
                     isgn = isgn + code(iorb)
                 enddo ! over iorb={1,norbs} loop
                 code(betta) = 0

             ! simulate one construct operator
             if (code(alpha) .eq. 0) then
                 do iorb=1,alpha-1
                     isgn = isgn + code(iorb)
                 enddo ! over iorb={1,norbs} loop
                 code(alpha) = 1

                 ! determine the row number and matrix element
                 dsgn = 1.0d0; if (mod(isgn, 2) .eq. 1) dsgn = - 1.0d0
                 jnew = jold - 2**(betta-1) + 2**(alpha-1); ibas = invsn(jnew)

                 if (ibas .eq. 0) stop "severe error happened in atomic_make_sp2np"
                 zmat_n(ibas, jbas) = zmat_n(ibas, jbas) + zmat_s(alpha, betta) * dsgn

             endif ! back if (code(betta) .eq. 0) block
             endif ! back if (code(alpha) .eq. 1) block

         enddo bettaloop ! over betta={1,norbs} loop
         enddo alphaloop ! over alpha={1,norbs} loop
     enddo jbasloop ! over jbas={1,ncfgs} loop

     return
  end subroutine atomic_make_sp2np

!-------------------------------------------------------------------------!
!>>> setup spin-orbit coupling matrix in {lz, sz} single particle basis <<!
!>>> norbs=02 => l=0; norbs=06 => l=1; norbs=10 => l=2; norbs=14 => l=3 <<!
!-------------------------------------------------------------------------!
  subroutine atomic_make_somatp( nband, nspin, norbs, lamb, somat)

     implicit none

! external arguments
! number of bands
     integer, intent(in) :: nband

! number of spins
     integer, intent(in) :: nspin

! number of orbits
     integer, intent(in) :: norbs

! spin-orbital coupling parameter
     real(8), intent(in) :: lamb

! spin-orbit coupling matrix in real-axis single particle basis
     complex(8), intent(out) :: somat(norbs, norbs)

! initialize somat to be zero
     somat = dcmplx(0.0D0, 0.0D0)

! spin orbit coupling in real orbit {px, sz}
     somat(1, 3) = - dcmplx(0.0D0, 1.0D0); somat(1, 6) = + dcmplx(1.0D0, 0.0D0)
     somat(2, 4) = + dcmplx(0.0D0, 1.0D0); somat(2, 5) = - dcmplx(1.0D0, 0.0D0)
     somat(3, 1) = + dcmplx(0.0D0, 1.0D0); somat(3, 6) = - dcmplx(0.0D0, 1.0D0)
     somat(4, 2) = - dcmplx(0.0D0, 1.0D0); somat(4, 5) = - dcmplx(0.0D0, 1.0D0)
     somat(5, 2) = - dcmplx(1.0D0, 0.0D0); somat(5, 4) = + dcmplx(0.0D0, 1.0D0)
     somat(6, 1) = + dcmplx(1.0D0, 0.0D0); somat(6, 3) = + dcmplx(0.0D0, 1.0D0)

     somat = somat * (lamb / 2.0D0)

     return
  end subroutine atomic_make_somatp
