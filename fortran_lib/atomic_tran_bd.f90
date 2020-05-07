  subroutine atomic_tran_blockdiag(ncfgs, c4mat_n, c4mat_bd, Umtrx, Umtrx1d, c4v_degs, num_2d)
    implicit none

! the number of configuration
    integer,      intent(in)  :: ncfgs

! the input target matrix     
    complex(8),   intent(in)  :: c4mat_n(ncfgs, ncfgs)

! the output block diagonal matrix    
    complex(8),   intent(out) :: c4mat_bd(ncfgs, ncfgs)

! the unitary matrix    
    complex(8),   intent(out) :: Umtrx(ncfgs, ncfgs)

! 1d version of unitary matrix
    integer,      intent(out) :: Umtrx1d(ncfgs)

! the output block diagonal matrix    
    integer   ,   intent(out) :: c4v_degs(ncfgs)

! the number of 2d representaiton
    integer   ,   intent(out) :: num_2d

! local variables
    integer     :: i, j, k, ibase, jbase, flag, tmp1, tmp2
    complex(8),   allocatable :: Umtrx1(:,:)
    complex(8),   allocatable :: Umtrx2(:,:)
    complex(8),   allocatable :: Umtrx3(:,:)
    complex(8),   allocatable :: c4mat_t(:,:)
    complex(8),   allocatable :: c4mat_t1(:,:)

! allocate
    allocate(Umtrx1(ncfgs, ncfgs))
    allocate(Umtrx2(ncfgs, ncfgs))
    allocate(Umtrx3(ncfgs, ncfgs))
    allocate(c4mat_t(ncfgs, ncfgs))
    allocate(c4mat_t1(ncfgs, ncfgs))

! initial
    ibase    = 0
    num_2d   = 0
    c4v_degs = 0
    Umtrx1   = dcmplx(0.d0, 0.d0)
    Umtrx2   = dcmplx(0.d0, 0.d0)
    Umtrx3   = dcmplx(0.d0, 0.d0)
    c4mat_t  = dcmplx(0.d0, 0.d0)
    c4mat_t1 = dcmplx(0.d0, 0.d0)
    do k = 1, ncfgs
        Umtrx1(k, k) = dcmplx(1.d0, 0.d0)
        Umtrx2(k, k) = dcmplx(1.d0, 0.d0)
        Umtrx3(k, k) = dcmplx(1.d0, 0.d0)
        Umtrx1d(k)   = k
    enddo
    c4mat_t = c4mat_n

! key part    
    i  = 1
    do while(i .le. ncfgs)
        ibase = 0
        Umtrx2   = dcmplx(0.d0, 0.d0)
        do k = 1, ncfgs
            Umtrx2(k, k) = dcmplx(1.d0, 0.d0)
        enddo
        do j = 1, ncfgs
            if(abs(c4mat_t(j,i)) .gt. 1E-6 .and. j .ge. i+2)then
                Umtrx2(i+1,  j) = dcmplx(1.d0, 0.d0)
                Umtrx2(j,  i+1) = dcmplx(1.d0, 0.d0)
                Umtrx2(i+1,i+1) = dcmplx(0.d0, 0.d0)
                Umtrx2(j,    j) = dcmplx(0.d0, 0.d0)
                ibase  = 2
                num_2d = num_2d +1
                tmp1         = Umtrx1d(i+1)
                Umtrx1d(i+1) = Umtrx1d(j)
                Umtrx1d(j)   = tmp1
                exit
            elseif(abs(c4mat_t(j,i)) .gt. 1E-6 .and. j .eq. i+1 )then
                ibase  = 2
                num_2d = num_2d +1
                exit
            elseif(abs(c4mat_t(j,i)) .gt. 1E-6 .and. i .eq. j )then
                ibase = 1
                exit
            endif
        enddo
        do j = i, i+ibase-1
            c4v_degs(j) = ibase - (j-i)
        enddo

        c4mat_t1 = matmul(c4mat_t, Umtrx2)
        c4mat_t  = matmul(transpose(Umtrx2),c4mat_t1)
        Umtrx1   = matmul(Umtrx1,Umtrx2)
        i = i + ibase
    enddo
    c4mat_bd = c4mat_t
    Umtrx    = Umtrx1
    do j = 1, ncfgs
        if(abs(Umtrx1(Umtrx1d(j),j)) .lt. 1e-6 )then
            call atomic_dump_isgnmat(ncfgs, Umtrx1d, 'atom_blockdiag_Umtrx1d.dat')
            call atomic_dump_c4mat2(ncfgs, Umtrx1,  'atom_blockdiag_Umtrx1.dat')
            print *,'j, Umtrx1d(j), Umtrx1(Umtrx1d(j),j)'
            print *, j, Umtrx1d(j), Umtrx1(Umtrx1d(j),j)
            stop "something wrong happens <atomic_tran_blockdiag-3>"
        endif
    enddo

! deallocate
    deallocate(Umtrx1)
    deallocate(Umtrx2)
    deallocate(Umtrx3)
    deallocate(c4mat_t)
    deallocate(c4mat_t1)
  end subroutine atomic_tran_blockdiag
!
!

  subroutine atomic_checkirrep(ncfgs, c4mat_bd, c2mat_bd, sgmvmat_bd, sgmdmat_bd, &
          c4v_degs, irrep_flag, irrep_type, num_rdu, num_irrep)
    implicit none

! the number of configuration
    integer,      intent(in)  :: ncfgs

! the output block diagonal matrix    
    complex(8),   intent(in)  :: c4mat_bd(ncfgs, ncfgs)

! the output block diagonal matrix    
    complex(8),   intent(in)  :: c2mat_bd(ncfgs, ncfgs)

! the output block diagonal matrix    
    complex(8),   intent(in)  :: sgmvmat_bd(ncfgs, ncfgs)

! the output block diagonal matrix    
    complex(8),   intent(in)  :: sgmdmat_bd(ncfgs, ncfgs)

! the output block diagonal matrix    
    integer   ,   intent(in)  :: c4v_degs(ncfgs)

! label every represenation is reduced or not   
    integer   ,   intent(out) :: irrep_flag(ncfgs)

! label every represenation is reduced or not   
    integer   ,   intent(out) :: irrep_type(ncfgs)

! the output block diagonal matrix    
    integer   ,   intent(out) :: num_rdu

! label the number every basis of every irrep.
! 1: A1, 2:A2, 3:B1, 4:B2, 5:E1_1, 6:E1_2
    integer   ,   intent(out) :: num_irrep(6)

! local variables
    integer     :: i, j, k, ibase, jbase, flag, num

! label whether reduced  or not
    integer     :: whether_irrep 

! mean the type of irrepresentation of C4v: 1 A1, 2 A2, 3 B1, 4 B2, 5 E    
    integer     :: irrep(2)

! contain the character of C4v in order, E,2C4,C2,Sigma_v,Sigma_d    
    real(8)     :: chrct(5)
    complex(8),   allocatable :: c4mat_t(:,:)
    complex(8),   allocatable :: c4mat_t1(:,:)

! allocate
    allocate(c4mat_t(ncfgs, ncfgs))
    allocate(c4mat_t1(ncfgs, ncfgs))

! initial
    irrep        = 0 
    irrep_flag   = 0
    irrep_type   = 0
    num_rdu      = 0
    ibase        = 0
    chrct        = 0.d0
    flag         = 0
    num_irrep    = 0

! key part    
    i = 1
    do while(i .le. ncfgs)
        num = c4v_degs(i)
        if (num .eq. 1)then
            chrct(1) = 1.d0
            chrct(2) = dble(c4mat_bd(i,i))
            chrct(3) = dble(c2mat_bd(i,i))
            chrct(4) = dble(sgmvmat_bd(i,i))
            chrct(5) = dble(sgmdmat_bd(i,i))
        elseif(num .eq. 2)then
            chrct(1) = 2.d0
            chrct(2) = dble(c4mat_bd(i,i)+c4mat_bd(i+1,i+1))
            chrct(3) = dble(c2mat_bd(i,i)+c2mat_bd(i+1,i+1))
            chrct(4) = dble(sgmvmat_bd(i,i) +sgmvmat_bd(i+1,i+1))
            chrct(5) = dble(sgmdmat_bd(i,i)+sgmdmat_bd(i+1,i+1))
        endif        
!        call atomic_whetherirrep(chrct, flag)
        call atomic_whetherirrep(chrct, whether_irrep, irrep)
!        write(71,*)chrct, whether_irrep, irrep
        do  j = i , i + num - 1
            irrep_flag(j) = whether_irrep
            irrep_type(j) = irrep(j- i+ 1)
            num_irrep(irrep(j- i+ 1)) = num_irrep(irrep(j- i+ 1)) + 1
        enddo
        i = i + num 
        if (.not. whether_irrep )then
            num_rdu = num_rdu + 1
        endif
    enddo
! average the two basis of E1 of C4v    
    num_irrep(5) = num_irrep(5) / 2
    num_irrep(6) = num_irrep(5)
    if(sum(num_irrep) .ne. ncfgs)then
        stop "the total number of basis of all irrep. goes wrong!<atomic_checkirrep>"
    endif

! deallocate
    deallocate(c4mat_t)
    deallocate(c4mat_t1)
  end subroutine atomic_checkirrep
!

  subroutine atomic_checkirrep2( ncfgs, c4v_umt, c4mat_bd, c2mat_bd, sgmvmat_bd, &
          sgmdmat_bd, num_irrep)
    implicit none

! the number of configuration
    integer,      intent(in)  :: ncfgs

! the transform matrix from Fock basis into irrep. of c4v 
    complex(8),   intent(in)  :: c4v_umt(ncfgs, ncfgs)

! the output block diagonal matrix    
    complex(8),   intent(in)  :: c4mat_bd(ncfgs, ncfgs)

! the output block diagonal matrix    
    complex(8),   intent(in)  :: c2mat_bd(ncfgs, ncfgs)

! the output block diagonal matrix    
    complex(8),   intent(in)  :: sgmvmat_bd(ncfgs, ncfgs)

! the output block diagonal matrix    
    complex(8),   intent(in)  :: sgmdmat_bd(ncfgs, ncfgs)

! label the number every basis of every irrep.
! 1: A1, 2:A2, 3:B1, 4:B2, 5:E1_1, 6:E1_2
    integer   ,   intent(in) :: num_irrep(6)

! local variables
    integer     :: i, j, k, ibase, jbase, flag, num

! label whether reduced  or not
    integer     :: whether_irrep 

! mean the type of irrepresentation of C4v: 1 A1, 2 A2, 3 B1, 4 B2, 5 E    
    integer     :: irrep(2)

! contain the character of C4v in order, E,2C4,C2,Sigma_v,Sigma_d    
    real(8)     :: chrct(5)
    complex(8),   allocatable :: c4mat_t(:,:)
    complex(8),   allocatable :: c4mat_t1(:,:)

! allocate
    allocate(c4mat_t(ncfgs, ncfgs))
    allocate(c4mat_t1(ncfgs, ncfgs))

! initial
    irrep        = 0 
    ibase        = 0
    chrct        = 0.d0
    flag         = 0

    num   =  sum(num_irrep(1: 4))
!    print *,'the total number of num_irrep of 1D is :', num 
! key part    
!    open(52,file = 'test_irrepright.dat', status = 'unknown')
!    write(52,'(A6,A10, A10,A10)')'No.','character','whether','irrep'
    i = 1
    do while(i .le. num)
        chrct(1) = 1.d0
        chrct(2) = dble(c4mat_bd(i,i))
        chrct(3) = dble(c2mat_bd(i,i))
        chrct(4) = dble(sgmvmat_bd(i,i))
        chrct(5) = dble(sgmdmat_bd(i,i))
!       call atomic_whetherirrep(chrct, flag)
        irrep = 0
        call atomic_whetherirrep(chrct, whether_irrep, irrep)
!        write(52,'(I5,2X,5I4,I5,I5)')i,chrct, whether_irrep, irrep(1)
        i = i + 1 
    enddo
!    close(52)
! average the two basis of E1 of C4v    

! deallocate
    deallocate(c4mat_t)
    deallocate(c4mat_t1)
  end subroutine atomic_checkirrep2

!
  subroutine atomic_whetherirrep(chrct, whether_irrep, irrep)
      implicit none
! the character of a certain representation of C4v
      real(8),      intent(in)  ::  chrct(5)

! whether reduced or not, 1: irreducible, 0: reducible
      integer,      intent(out) ::  whether_irrep

! the kind of irrep of C4v: 1 A1, 2 A2, 3 B1, 4 B2, 5 E
      integer,      intent(out) ::  irrep(2)

! local variables 
      integer    :: i, j, k
      real(8)    :: CHA(5,5), weight(5)
      real(8)    :: cht_t, cht, tmp1, tmp2
      integer    :: cnt


! initial
      CHA           =  0.d0
      cht_t         =  0.d0
      cht           =  0.d0
      tmp1          =  1.d0
      tmp2          =  1.d0
      whether_irrep =  5.d0
      irrep         =  0

! key part
!     define the character table of C4v 
      CHA(1,:)  = (/ 1.d0, 1.d0, 1.d0, 1.d0, 1.d0 /)
      CHA(2,:)  = (/ 1.d0, 1.d0, 1.d0,-1.d0,-1.d0 /)
      CHA(3,:)  = (/ 1.d0,-1.d0, 1.d0, 1.d0,-1.d0 /)
      CHA(4,:)  = (/ 1.d0,-1.d0, 1.d0,-1.d0, 1.d0 /)
      CHA(5,:)  = (/ 2.d0, 0.d0,-2.d0, 0.d0, 0.d0 /)

      weight(:) = (/ 1.d0, 2.d0, 1.d0, 2.d0, 2.d0 /) 
      cht_t = 0.d0
      do j = 1, 5
          cht_t = cht_t + weight(j)*chrct(j)**2
      enddo
      if (abs(cht_t - 8.d0) .lt. 1e-6 )then
          whether_irrep = 1
      elseif(abs(cht_t - 16.d0) .lt. 1E-6)then
          whether_irrep = 0
      else
          print *, "cht_t", cht_t
          write(70,*) "chrct",chrct
          stop "something wrong happens in atomic_whetherirrep<1>"
      endif
      if(whether_irrep .eq. 1)then
          if(abs(chrct(1)- 2.d0) .lt. 1e-6)then
              irrep(1) = 5
              irrep(2) = 5
          elseif(abs(chrct(1) - 1.d0) .lt. 1e-6)then
              do k =1, 4
                  tmp1 = 2.d0
                  call rdot_product(5, weight, chrct, CHA(k,:), tmp1) 
                  if(abs(tmp1 -  8.d0) .lt. 1E-6 .or. abs(tmp1 + 8.d0) .lt. 1e-6)then
                      irrep(1) = k
                      exit
                  elseif(abs(tmp1) .lt. 1e-6 )then
                      continue
                  elseif(abs(tmp1) .lt. 1e-6 .and. k .eq. 4)then
                      stop "something happen in atomic_whetherirrep<4>"
                  else
                      stop "something happen in atomic_whetherirrep<3>"
                  endif
              enddo
          else
              stop "something wrong happen in atomic_whetherirrep<2>"
          endif
      elseif(whether_irrep .eq. 0)then
! reducing           
          if(abs(chrct(1) - 2.d0) .gt. 1e-6)then
              stop "something wrong happens in atomic_whetherirrep<6> because &
              1D rep. can't be reduced"
          else
              cnt = 0
              do k = 1, 4
                  tmp2 = 5.d0
                  call rdot_product(5, weight, chrct, CHA(k, :), tmp2)
                  if(abs(tmp2 - 8.d0) .lt. 1e-6)then
                      cnt = cnt +1
                      irrep(cnt) = k
                  endif
                  if(cnt .eq. 2)then
                      exit 
                  endif
              enddo
          endif
      else
          stop "something wrong happen in atomic_whetherirrep<5>"
      endif

  end subroutine atomic_whetherirrep
! 
  subroutine rdot_product(num,weight,a,b,c)
      implicit none
      integer,      intent(in)  :: num
      real(8),      intent(in)  :: weight(num)
      real(8),      intent(in)  :: a(num)
      real(8),      intent(in)  :: b(num)
      real(8),      intent(out) :: c

!
      integer  :: i
      c = 0.0
      do i = 1, num
          c = c + a(i)*b(i)*weight(i)
      enddo

      return
  end subroutine rdot_product
!
!
  subroutine atomic_take_reduce(ncfgs, norbs, c4v_umt_1d, irrep_flag, irrep_type, &
          invcd, basis, c4mat_bd, c2mat_bd, sgmvmat_bd, sgmdmat_bd, num_rdu)
      implicit none
      
! the number of configurations      
      integer,       intent(in)   ::  ncfgs

! the number of configurations      
      integer,       intent(in)   ::  norbs

! the 1d version of the unitary matrix of c4v group
      integer,       intent(in)   ::  c4v_umt_1d(ncfgs)

! label every represenation is reduced or not   
      integer,       intent(in)   :: irrep_flag(ncfgs)

! label every represenation's kind
      integer,       intent(in)   :: irrep_type(ncfgs)

! the output block diagonal matrix    
      complex(8),    intent(in)   :: c4mat_bd(ncfgs, ncfgs)

! the decimal number of Fock basis
      integer,       intent(in)   :: basis(ncfgs)

! the binary number of Fock basis
      integer,       intent(in)   :: invcd(norbs,ncfgs)

! the output block diagonal matrix    
      complex(8),    intent(in)   :: c2mat_bd(ncfgs, ncfgs)

! the output block diagonal matrix    
      complex(8),    intent(in)   :: sgmvmat_bd(ncfgs, ncfgs)

! the output block diagonal matrix    
      complex(8),    intent(in)   :: sgmdmat_bd(ncfgs, ncfgs)

! the number of configurations      
      integer,       intent(in)   ::  num_rdu

! local variables 
      integer    ::   i, j, k, cnt, ibase, jbase, cnt2

! initial
      ibase  = 0
      cnt    = 0
      cnt2   = 0
! key part
      i = 1
      open(44, file = "atom.reducible.dat", status = "unknown")
      write(44,'(A10,5x,A15,4x,A10,4x,A10,4x,A10,4x,A10)') 'Number','Fock state',&
          'C4', 'C2','sigma_v','sigma_d'
      do while( i .le. ncfgs)
          ibase = 0
          if(irrep_flag(i) .eq. 0)then
              ibase = 2
              cnt   = cnt + 1 
              write(44,'(I3,5X,10I1,4x,2f5.1,4x,2f5.1,4x,2f5.1,4x,2f5.1)') cnt, invcd(1:norbs,c4v_umt_1d(i)), &
                  dble(c4mat_bd(i,i)),dble(c4mat_bd(i,i+1)),dble(c2mat_bd(i,i)),dble(c2mat_bd(i,i+1)), &
                  dble(sgmvmat_bd(i,i)),dble(sgmvmat_bd(i,i+1)),dble(sgmdmat_bd(i,i)),dble(sgmdmat_bd(i,i+1))
              write(44,'(I3,5X,10I1,4x,2f5.1,4x,2f5.1,4x,2f5.1,4x,2f5.1)') cnt, invcd(1:norbs,c4v_umt_1d(i+1)), &
                  dble(c4mat_bd(i+1,i)),dble(c4mat_bd(i+1,i+1)),dble(c2mat_bd(i+1,i)),dble(c2mat_bd(i+1,i+1)), &
                  dble(sgmvmat_bd(i+1,i)),dble(sgmvmat_bd(i+1,i+1)),dble(sgmdmat_bd(i+1,i)),dble(sgmdmat_bd(i+1,i+1))
          else
              ibase = 1
          endif
          i = i + ibase
      enddo
      close(44)

      if(cnt .ne. num_rdu)then
          stop "the number of reducible is not consistent with the input data"
      endif

! for other purpose      
      i = 1
      open(44, file = "atom.irrep_matrix.dat", status = "unknown")
      write(44,'(A10,5x,A15,4x,A10,4x,A10,4x,A10,4x,A10)') 'Number','Fock state',&
          'C4', 'C2','sigma_v','sigma_d'
      do while( i .le. ncfgs)
          ibase = 0
          if(irrep_flag(i) .eq. 1 .and. irrep_type(i) .eq. 5)then
              ibase = 2
              cnt2   = cnt2 + 1 
              write(44,'(5x,I3,5X,10I1,4x,2f5.1,4x,2f5.1,4x,2f5.1,4x,2f5.1)') i, invcd(1:norbs,c4v_umt_1d(i)), &
                  dble(c4mat_bd(i,i)),dble(c4mat_bd(i,i+1)),dble(c2mat_bd(i,i)),dble(c2mat_bd(i,i+1)), &
                  dble(sgmvmat_bd(i,i)),dble(sgmvmat_bd(i,i+1)),dble(sgmdmat_bd(i,i)),dble(sgmdmat_bd(i,i+1))
              write(44,'(5x,I3,5X,10I1,4x,2f5.1,4x,2f5.1,4x,2f5.1,4x,2f5.1)') i+1, invcd(1:norbs,c4v_umt_1d(i+1)), &
                  dble(c4mat_bd(i+1,i)),dble(c4mat_bd(i+1,i+1)),dble(c2mat_bd(i+1,i)),dble(c2mat_bd(i+1,i+1)), &
                  dble(sgmvmat_bd(i+1,i)),dble(sgmvmat_bd(i+1,i+1)),dble(sgmdmat_bd(i+1,i)),dble(sgmdmat_bd(i+1,i+1))
          else
              ibase = 1
          endif
          i = i + ibase
      enddo
      close(44)
! other purpose end      
  end subroutine atomic_take_reduce
!
  subroutine atomic_reduction_unitary(A)
      implicit none
      complex(8),  intent(out)  :: A(2,2)

! local variables
      real(8)     :: sqrt2, sqrt2m
      complex(8)  :: cone, czero, cima
! definition
      cone     = dcmplx(1.d0,  0.d0)
      czero    = dcmplx(0.d0,  0.d0)
      cima     = dcmplx(0.d0,  1.d0)
      sqrt2m   = 1.d0/sqrt(2.d0) 

! key part      
      A(1, :) = sqrt2m*(/ cone ,  cone  /)
      A(2, :) = sqrt2m*(/ cone , -cone  /)
  end subroutine atomic_reduction_unitary
!
  subroutine atomic_irrep_umtrx(ncfgs, irrep_flag, irrep_type, num_irrep, c4degs,& 
      rdct_umt, c4mat_bd, c2mat_bd, sgmvmat_bd, sgmdmat_bd, Umtrx,  c4virrepdegs)
      implicit none
      
! the number of configurations      
      integer,       intent(in)     ::  ncfgs

! number of basis of every irrep.      
      integer,       intent(in)     ::  num_irrep(6)

! label every represenation is reduced or not   
      integer,       intent(in)     ::  c4degs(ncfgs)

! label every represenation is reduced or not   
      integer,       intent(inout)  ::  irrep_flag(ncfgs)

! label every represenation is reduced or not   
      integer,       intent(inout)  ::  irrep_type(ncfgs)

! reduce the 2d reducible block matrix  
      complex(8),    intent(in)     ::  rdct_umt(2,2)

! the output block diagonal matrix    
      complex(8),    intent(in)     ::  c4mat_bd(ncfgs, ncfgs)

! the output block diagonal matrix    
      complex(8),    intent(in)     ::  c2mat_bd(ncfgs, ncfgs)

! the output block diagonal matrix    
      complex(8),    intent(in)     ::  sgmvmat_bd(ncfgs, ncfgs)

! the output block diagonal matrix    
      complex(8),    intent(in)     ::  sgmdmat_bd(ncfgs, ncfgs)

! the unitary matrix    
      complex(8),    intent(inout)  ::  Umtrx(ncfgs, ncfgs)


! the output block diagonal matrix    
      integer   ,    intent(out)    ::   c4virrepdegs(ncfgs)

! local variables
! the output block diagonal matrix    
      integer,     allocatable :: c4degs1(:)
! label every represenation is reduced or not   
      integer,     allocatable :: irrep_flag1(:)
! label every represenation is reduced or not   
      integer,     allocatable :: irrep_type1(:)
! the output block diagonal matrix    
      integer                  :: num_rdu
! label the number every basis of every irrep.
! 1: A1, 2:A2, 3:B1, 4:B2, 5:E1_1, 6:E1_2
      integer                  :: num_irrep1(6)

      integer     ::  i, j, k, ibase,ibase1(6), ibase2(6)
      real(8)     ::  tmp1, tmp2
      complex(8)  :: cone, czero, cima
      complex(8), allocatable  ::  Umtrx1(:, :), Umtrx2(:, :)
      complex(8), allocatable  ::  Umtrx3(:, :), Umtrx4(:, :)
      complex(8), allocatable  ::  c4mat_bd2(:, :)
      complex(8), allocatable  ::  c4mat_irrep(:, :)
      complex(8), allocatable  ::  c2mat_irrep(:, :)
      complex(8), allocatable  ::  sgmvmat_irrep(:, :)
      complex(8), allocatable  ::  sgmdmat_irrep(:, :)

! allocate
      allocate(c4degs1(ncfgs))
      allocate(irrep_flag1(ncfgs))
      allocate(irrep_type1(ncfgs))
      allocate(Umtrx1(ncfgs,ncfgs))
      allocate(Umtrx2(ncfgs,ncfgs))
      allocate(Umtrx3(ncfgs,ncfgs))
      allocate(Umtrx4(ncfgs,ncfgs))
      allocate(c4mat_bd2(ncfgs,ncfgs))
      allocate(c4mat_irrep(ncfgs,ncfgs))
      allocate(c2mat_irrep(ncfgs,ncfgs))
      allocate(sgmvmat_irrep(ncfgs,ncfgs))
      allocate(sgmdmat_irrep(ncfgs,ncfgs))

! definition
      cone     = dcmplx(1.d0,  0.d0)
      czero    = dcmplx(0.d0,  0.d0)
      cima     = dcmplx(0.d0,  1.d0)

      
! initial
      Umtrx1       = czero
      Umtrx2       = czero
      Umtrx3       = czero
      Umtrx4       = czero
      ibase        = 0
      irrep_flag1  = 0
      irrep_type1  = 0
      num_irrep1   = 0
      ibase2       = 0
      ibase1(1)    = 0
      ibase1(2)    = ibase1(1)+ num_irrep(1)
      ibase1(3)    = ibase1(2)+ num_irrep(2)
      ibase1(4)    = ibase1(3)+ num_irrep(3)
      ibase1(5)    = ibase1(4)+ num_irrep(4)
      ibase1(6)    = ibase1(5)+ num_irrep(5)
      ibase2       = ibase1 + ibase2
      c4virrepdegs = 0
      c4mat_bd2    = c4mat_bd
      c4degs1      = c4degs

      Umtrx3       = Umtrx
      Umtrx4       = Umtrx
! key part
      do i = 1, ncfgs
          Umtrx1( i, i) = cone
      enddo
 
      j = 1 
      do while(j .le. ncfgs)
          ibase = 0
          if(irrep_flag(j) .eq. 0)then
              Umtrx1(j,     j) = rdct_umt(1,1)
              Umtrx1(j,   j+1) = rdct_umt(1,2)
              Umtrx1(j+1,   j) = rdct_umt(2,1)
              Umtrx1(j+1, j+1) = rdct_umt(2,2)
              c4degs1(j)       = 1
              c4degs1(j+1)     = 1
              ibase = 2
!              irrep_flag1(j)   = 1
!              irrep_flag1(j+1) = 1
          elseif(irrep_flag(j) .eq. 1)then
              ibase = 1
!              irrep_flag1(j)   = 1
          else
              stop "something wrong happens in <atomic_irrep_umtrx>"
          endif
          j = j + ibase
      enddo
      Umtrx4 = matmul(Umtrx3, Umtrx1)
      call atomic_dump_c4mat2(ncfgs, Umtrx1, 'atom.bd-irrep.dat')

      call atomic_tran_c4mat( ncfgs, Umtrx1,  c4mat_bd,   c4mat_irrep) 
      call atomic_tran_c4mat( ncfgs, Umtrx1,  c2mat_bd,   c2mat_irrep) 
      call atomic_tran_c4mat( ncfgs, Umtrx1,  sgmvmat_bd, sgmvmat_irrep) 
      call atomic_tran_c4mat( ncfgs, Umtrx1,  sgmdmat_bd, sgmdmat_irrep) 

      call atomic_checkirrep(ncfgs, c4mat_irrep, c2mat_irrep, sgmvmat_irrep, &
      sgmdmat_irrep, c4degs1, irrep_flag1, irrep_type1, num_rdu, num_irrep1)

!      
      j = 1
      do while(j .le. ncfgs)
          ibase = 0
          if(irrep_type1(j) .eq. 5)then
              if(abs(dble(c4mat_bd2(j+1,j)) - 1.d0) .lt. 1e-6)then
                  ibase1(5) = ibase1(5) + 1
                  ibase1(6) = ibase1(6) + 1
                  Umtrx2(j,ibase1(5)) = cone
                  Umtrx2(j+1,ibase1(6)) = cone
                  ibase =  2
              elseif(abs(dble(c4mat_bd2(j+1,j)) + 1.d0) .lt. 1e-6)then
                  ibase1(5) = ibase1(5) + 1
                  ibase1(6) = ibase1(6) + 1
                  Umtrx2(j+1,ibase1(5)) = cone
                  Umtrx2(j,ibase1(6)) = cone
                  ibase =  2
              endif
          else
              do i = 1, 4
                  if(irrep_type1(j) .eq. i )then
                      ibase1(i) = ibase1(i) + 1
                      Umtrx2(j, ibase1(i)) = cone
                      ibase = 1
                      exit
                  endif
              enddo
          endif
          j = j+ ibase
      enddo
      call atomic_dump_c4mat2(ncfgs, Umtrx2, 'atom.irrep-degs.dat')
      Umtrx = matmul(Umtrx4, Umtrx2)

      irrep_flag = irrep_flag1


      k = 0
      do k = 1, 6
          if(k .lt. 6)then
              do i = 1, num_irrep(k)
                  c4virrepdegs(ibase2(k+1)-i+1) = i
                  irrep_type(ibase2(k)+i)       = k
              enddo
          elseif(k .eq. 6)then
              do i = 1, num_irrep(k)
                  c4virrepdegs(ibase2(k)+num_irrep(6)-i+1) = i
                  irrep_type(ibase2(k)+i)       = 5
              enddo
          else
              stop "something wrong happen in atomic_irrep_umtrx-3"
          endif
      enddo

! deallocate
      deallocate(Umtrx1)
      deallocate(Umtrx2)
      deallocate(Umtrx3)
      deallocate(Umtrx4)
      deallocate(c4mat_bd2)
      deallocate(c4mat_irrep)
      deallocate(c2mat_irrep)
      deallocate(sgmvmat_irrep)
      deallocate(sgmdmat_irrep)
      deallocate(c4degs1)
      deallocate(irrep_flag1)
      deallocate(irrep_type1)
  end subroutine atomic_irrep_umtrx
! 
! aim : to calculate the number of vpms based on a degeneracy vec
  subroutine count_vpm(ncfgs, deg_vec, cnt, list)
      implicit none
      integer,      intent(in)   ::  ncfgs

! the degeneracy information      
      integer,      intent(in)   ::  deg_vec(ncfgs)

! the number of vpms      
      integer,      intent(out)  ::  cnt
      
! the number of vpms      
      integer,      intent(out)  ::  list(ncfgs)
      
! local variables
      integer   ::  i, ibase

! initial 
      ibase    = 0
      cnt      = 0
      list     = 0
      i = 1
      do while(i .le. ncfgs)
          cnt = cnt + deg_vec(i)**2
          list(deg_vec(i)) = list(deg_vec(i)) + 1
          ibase = deg_vec(i)
          i = i + ibase
      end do

  end subroutine count_vpm

