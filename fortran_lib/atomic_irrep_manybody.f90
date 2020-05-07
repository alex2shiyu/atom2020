!
!
!
  subroutine make_irrep_mb(norbs,ncfgs,ntots,nstat,ntiny,nlarge,Dn_num,Dn_dim,Dn_vec)
      implicit none

! the number of orbitals(spin included)      
      integer,   intent(in)      :: norbs 

! the number of configurations      
      integer,   intent(in)      :: ncfgs 

! the lower boundry of paritcle number
      integer,   intent(in)      :: ntots

! the lower boundry of paritcle number
      integer,   intent(in)      :: nstat(0:norbs) 

! the lower boundry of paritcle number
      integer,   intent(in)      :: ntiny 

! the top boundry of paritcle number
      integer,   intent(in)      :: nlarge

! the number of set of basis of irreps of each type, for example in the N=3 subspace, there
! maybe many pairs of many body basis which can span this type of irrep
      integer,   intent(out)     :: Dn_num(5)

! store the dimensionality of every irreps
      integer,   intent(out)     :: Dn_dim(5)

! store the key information: number of type of each irrep; the number of set of
! basis of irrep of each type(D1_num); label of each basis of irrep, maximally 2
! of C4v; the basis vector of each basis of each set of each irrep in the Fock
! representation
      complex(8),   intent(out)     :: Dn_vec(5,1000,2,ncfgs)

! temparay variables
! the decimal number(2^0, ...) of each states(1, 2, 3, 4, ...)
      integer     :: basis(ncfgs)

! the serials number of each decimal number, contrary to <basis>
      integer     :: invsn(0:2**norbs-1)

! the binary information of each serial number state
      integer     :: invcd(norbs, ncfgs)

! maybe many pairs of many body basis which can span this type of irrep
      integer     :: D1_num(5)

! store the dimensionality of every irreps
      integer     :: D1_dim(5)

! store the key information: number of type of each irrep; the number of set of
! basis of irrep of each type(D1_num); label of each basis of irrep, maximally 2
! of C4v; the basis vector of each basis of each set of each irrep in the Fock
! representation
      complex(8),   intent(out)     :: D1_vec(5,1000,2,ncfgs)



! The auxiliary matrix
      integer     :: D1_code(5,1000,2,norbs)

! local variables
      integer     :: i, j, k, ii, jj, kk
      integer     :: basis1, basis2, serial1(norbs)
      

! construct Fock basis sheet (just follow nhtong's note)
      call atomic_make_basis( norbs, ncfgs, ntots, nstat, basis, invcd, invsn)
      call make_D1(norbs,ncfgs,basis,invsn,invcd,D1_num,D1_dim,D1_vec)


  end subroutine make_irrep_mb

  subroutine make_D1(norbs,ncfgs,basis,invsn,invcd,D1_num,D1_dim,D1_vec)
      implicit none

! the number of orbitals(spin included)      
      integer,   intent(in)      :: norbs 

! the number of configurations      
      integer,   intent(in)      :: ncfgs 

! the decimal number(2^0, ...) of each states(1, 2, 3, 4, ...)
      integer,   intent(in)      :: basis(ncfgs)

! the serials number of each decimal number, contrary to <basis>
      integer,   intent(in)      :: invsn(0:2**norbs-1)

! the binary information of each serial number state
      integer,   intent(in)      :: invcd(norbs, ncfgs)

! the number of set of basis of irreps of each type, for example in the N=3 subspace, there
! maybe many pairs of many body basis which can span this type of irrep
      integer,   intent(out)     :: D1_num(5)

! store the dimensionality of every irreps
      integer,   intent(out)     :: D1_dim(5)

! store the key information: number of type of each irrep; the number of set of
! basis of irrep of each type(D1_num); label of each basis of irrep, maximally 2
! of C4v; the basis vector of each basis of each set of each irrep in the Fock
! representation
      complex(8),   intent(out)     :: D1_vec(5,1000,2,ncfgs)

! The auxiliary matrix
      integer  :: D1_code(5,1000,2,norbs)

! local variables
      integer     :: i, j, k, ii, jj, kk
      integer     :: basis1, basis2, code1(norbs), code2(norbs)
      integer     :: dnum1,dnum2

! initial 
      basis2  = 0
      D1_code = 0
      D1_vec  = 0
      code1   = 0
      code2   = 0

! set D1_num      
      do i = 1, 4
          D1_num(i) = 2
          D1_dim(i) = 1
      enddo
      D1_num(5) = 2
      D1_dim(5) = 2
      D1_code(1, 1, 1, 7) = 1
      D1_code(1, 2, 1, 8) = 1
      D1_code(3, 1, 1, 9) = 1
      D1_code(3, 2, 1,10) = 1
      D1_code(4, 1, 1, 5) = 1
      D1_code(4, 2, 1, 6) = 1
      D1_code(5, 1, 1, 1) = 1
      D1_code(5, 1, 2, 3) = 1
      D1_code(5, 2, 1, 2) = 1
      D1_code(5, 2, 2, 4) = 1

      do i = 1, 5
          do j = 1, D1_num(i)
              do k = 1, D1_dim(i)
                  do ii = 1, norbs
                      basis2 = basis2 + 2**(D1_code(i,j,k,ii))
                  enddo
                  D1_vec(i,j,k, invsn(basis2)) = dcmplx(1.d0,0.d0) 
              enddo
          enddo
      enddo
  end subroutine make_D1

  subroutine make_Dn(norbs,ncfgs,basis,invsn,invcd,D2_num,D2_dim,D2_vec,D1_num,D1_dim,D1_vec,Dn_num,Dn_dim,Dn_vec)
      implicit none

! the number of orbitals(spin included)      
      integer,   intent(in)      :: norbs 

! the number of configurations      
      integer,   intent(in)      :: ncfgs 

! the decimal number(2^0, ...) of each states(1, 2, 3, 4, ...)
      integer,   intent(in)      :: basis(ncfgs)

! the serials number of each decimal number, contrary to <basis>
      integer,   intent(in)      :: invsn(0:2**norbs-1)

! the binary information of each serial number state
      integer,   intent(in)      :: invcd(norbs, ncfgs)

! the number of set of basis of irreps of each type, for example in the N=3 subspace, there
! maybe many pairs of many body basis which can span this type of irrep
      integer,   intent(in)      :: D1_num(5)

! store the dimensionality of every irreps
      integer,   intent(in)      :: D1_dim(5)

! store the key information: number of type of each irrep; the number of set of
! basis of irrep of each type(D1_num); label of each basis of irrep, maximally 2
! of C4v; the basis vector of each basis of each set of each irrep in the Fock
! representation
      complex(8),   intent(in)      :: D1_vec(5,1000,2,ncfgs)

! the number of set of basis of irreps of each type, for example in the N=3 subspace, there
! maybe many pairs of many body basis which can span this type of irrep
      integer,   intent(in)      :: D2_num(5)

! store the dimensionality of every irreps
      integer,   intent(in)      :: D2_dim(5)

! store the key information: number of type of each irrep; the number of set of
! basis of irrep of each type(D1_num); label of each basis of irrep, maximally 2
! of C4v; the basis vector of each basis of each set of each irrep in the Fock
! representation
      complex(8),   intent(in)      :: D2_vec(5,1000,2,ncfgs)

! the number of set of basis of irreps of each type, for example in the N=3 subspace, there
! maybe many pairs of many body basis which can span this type of irrep
      integer,   intent(out)     :: Dn_num(5)

! store the dimensionality of every irreps
      integer,   intent(out)     :: Dn_dim(5)

! store the key information: number of type of each irrep; the number of set of
! basis of irrep of each type(D1_num); label of each basis of irrep, maximally 2
! of C4v; the basis vector of each basis of each set of each irrep in the Fock
! representation
      complex(8),   intent(out)     :: Dn_vec(5,1000,2,ncfgs)

! this part will prepare the CG efficient of two irrep
! the number of reduction sub irreps when ith irrep and jth irrep direct product
      integer     :: dp_num(5,5)

! the dimension of every reduction irrep , the number of reduction sub-irrep is 4 at most
      integer     :: dp_type(5,5,4)

! the Cg Coefficients: the number of one direct product matrix; another direct
! product matrix; the number of reduction matrix; the row label; the column
! labels
      complex(8)  :: dp_vec(5,5,4,4,2)

! The auxiliary matrix
      integer     :: D1_code(5,1000,2,norbs)
      complex(8)  :: dpterm_vec(ncfgs,4)

! local variables
      integer     :: i, j, k, ii, jj, kk, n, nn, i1, i2, i3
      integer     :: basis1, basis2, serial1(norbs),code1(norbs)
      integer     :: code2(norbs),code3(norbs),coden(norbs,100),code_new(norbs)
      integer     :: sign1,sign2,sign3,sign4,dnum1,dnum2,flag1,flag2,flag3
      integer     :: sn1,sn2,newdnm
      integer     :: vec1(ncfgs),vec2(ncfgs),vec3(ncfgs),vec4(ncfgs)
      integer     :: t1,t2,cnt1,cnt2
      complex(8)  :: tmp1, tmp2

! initial
      flag2      = 0
      vec1       = 0
      vec2       = 0
      coden      = 0
      cnt1       = 0
      cnt2       = 0
      code_new   = 0
      dpterm_vec = 0
! key part       
      
! make Cg coefficient
      call make_Cg(D1_dim,dp_num,dp_type,dp_vec)
      do i = 1, 5
          do j = 1, D2_num(i)
              do ii = 1, 5
                  do jj = 1, D1_num(ii)
!                     
! check whether these two irrepresentation can direct product,
! which  means c^\dag_2|0100> is not allowed
                      !-------------------
                      call check_directproduct(norbs,ncfgs,basis,invsn,invcd,D1_dim,&
                          D1_vec(ii,jj,:,:),D2_dim,D2_vec(i,j,:,:),flag2)
! check end
                      ! if direct product is valid, this part will work which is
                      ! the key part of this "atomic_irrep_manybody.f90"
!-------------------------when direct product is possible,define joint basis,now-------------                      
                      if(flag2 .eq. 0)then
                          cnt1 = 0
                          dpterm_vec = 0
                          do k = 1, D2_dim(i)
                              do kk = 1, D1_dim(ii)
                                  cnt1 = cnt1 + 1
                                  vec3 = 0
                                  do n = 1, ncfgs
                                      code3 = 0
                                      if (D2_vec(i,j,k,n) .ne. 0)then
                                          tmp1 = D2_vec(i,j,k,n)
                                          t2 = invsn(n)
                                          code3 = invcd(t2)
                                          !----------------
                                          do nn = 1, ncfgs
                                              if( D1_vec(ii, jj, kk, nn) .ne.0)then
                                                  tmp2 = D1_vec(ii, jj, kk, nn)
                                                  t3   = invsn(nn) 
                                                  call get_power(nn, 2, 10, sn2)
                                                  exit
                                              endif
                                          enddo
                                          !----------------
                                          call get_newstate(norbs,ncfgs,basis,invsn,invcd,code3,sn2,sign4,code_new,newdnm)
                                          vec3(newdnm) = sign4*tmp2*tmp1
                                      endif
                                  enddo
                                  dpterm_vec(:,cnt1) = vec3
!-------------------------joint basis defined --------------------------------                                  
!
                              enddo
                          enddo ! loop over dimension of a certain irrep of D2!
                      endif
!-----------------------------------define new basis in whole space using CG coeficients--------------------------- 
                      if(cnt1 .ne. dp_num(i,j))
!
                  enddo
              enddo
          enddo
      enddo

  end subroutine make_Dn

  subroutine check_directproduct(norbs,ncfgs,basis,invsn,invcd,D1_dim,D1_vec,D2_dim,D2_vec,flag)
      implicit none
! the number of orbitals(spin included)      
      integer,   intent(in)      :: norbs 

! the number of configurations      
      integer,   intent(in)      :: ncfgs 

! the decimal number(2^0, ...) of each states(1, 2, 3, 4, ...)
      integer,   intent(in)      :: basis(ncfgs)

! the serials number of each decimal number, contrary to <basis>
      integer,   intent(in)      :: invsn(0:2**norbs-1)

! the binary information of each serial number state
      integer,   intent(in)      :: invcd(norbs, ncfgs)

! store the dimensionality of every irreps
      integer,   intent(in)      :: D1_dim

! store the key information: number of type of each irrep; the number of set of
! basis of irrep of each type(D1_num); label of each basis of irrep, maximally 2
! of C4v; the basis vector of each basis of each set of each irrep in the Fock
! representation
      complex(8),   intent(in)   :: D1_vec(2,ncfgs)

! store the dimensionality of every irreps
      integer,   intent(in)      :: D2_dim

! store the key information: number of type of each irrep; the number of set of
! basis of irrep of each type(D1_num); label of each basis of irrep, maximally 2
! of C4v; the basis vector of each basis of each set of each irrep in the Fock
! representation
      complex(8),   intent(in)   :: D2_vec(2,ncfgs)

      integer,   intent(out)     :: flag

! the Cg Coefficients: the number of one direct product matrix; another direct
! product matrix; the number of reduction matrix; the row label; the column
! labels
! The auxiliary matrix
      integer     :: D1_code(5,1000,2,norbs)
      complex(8)  :: dpterm_vec(ncfgs,4)

! local variables
      integer     :: i, j, k, ii, jj, kk, n, nn, i1, i2, i3
      integer     :: basis1, basis2, serial1(norbs),code1(norbs)
      integer     :: code2(norbs),code3(norbs),coden(norbs,100),code_new(norbs)
      integer     :: sign1,sign2,sign3,sign4,dnum1,dnum2,flag1,flag2,flag3
      integer     :: sn1,sn2,newdnm
      integer     :: vec1(ncfgs),vec2(ncfgs),vec3(ncfgs),vec4(ncfgs)
      integer     :: t1,t2,cnt1,cnt2
      complex(8)  :: tmp1, tmp2

! initial
      vec1       = 0
      vec2       = 0
      coden      = 0
      cnt1       = 0
      cnt2       = 0
      code_new   = 0
      dpterm_vec = 0
      code1      = 0
! key part       
      do k = 1, D2_dim
          do n = 1, ncfgs
              if( abs(D2_vec(k,n)) .gt. 1E-6)then
                  t1      = invsn(n)
                  serial1 = invcd(t1)
                  do nn = 1, norbs
                      if(serial1(nn) .eq. 1)then
                          code1(nn) = 1
                      endif
                  enddo
              endif
          enddo
      enddo
      !---------------------
      flag2 = 0
      do kk = 1, D1_dim
          dnum1 = 0
          do nn = 1, ncfgs
              if(abs(D1_vec(kk,nn)) .gt. 1E-7)then
                  dnum1 = nn
                  exit
              endif
          enddo
          flag1 = 0
          do nn = 1, norbs
              if(2**(nn-1) .eq. dnum1)then
                  if( code1(nn) .eq. 1)then
                      flag2 = 1
                      exit
                  else
                      exit
                  endif
              elseif(nn .eq. norbs)then
                  stop "something wrong happens here in make_Dn!"
              else
                  continue
              endif
          enddo
          if(flag2 .eq. 1)then
              exit
          endif
      enddo ! end loop D1_dim
      flag = flag2

      return
  end subroutine check_directproduct


  subroutine get_newstate(norbs,ncfgs,basis,invsn,invcd,code,sn,sgn_t,code_new,dnum)
      implicit none

! the number of orbitals(spin included)      
      integer,   intent(in)      :: norbs 

! the number of configurations      
      integer,   intent(in)      :: ncfgs 

! the decimal number(2^0, ...) of each states(1, 2, 3, 4, ...)
      integer,   intent(in)      :: basis(ncfgs)

! the serials number of each decimal number, contrary to <basis>
      integer,   intent(in)      :: invsn(0:2**norbs-1)

! the binary information of each serial number state
      integer,   intent(in)      :: invcd(norbs, ncfgs)

! the binary code for old states
      integer,   intent(in)      :: code(norbs)
  
! the binary code for old states
      integer,   intent(in)      :: sn
  
! the binary code for old states
      integer,   intent(out)      :: code_new(norbs)

! the sign of new states due to anti-commute relation
      integer,   intent(out)      :: sgn_t
 
! the value of binary type of new state
      integer,   intent(out)      :: dnum

! local variables 
      integer      ::   i, j, k, sgn
      integer      ::   code_t(norbs)

! key part

! first check 
      if(code(sn) .eq. 1)then
          stop "wrong happens in get_newstate!"
      endif
      sgn  = 1
      do i = 1, sn-1
          if(code(i) .eq. 1)then
              sgn = sgn * (-1)
          endif
      enddo
      sgn_t    = sgn
      code(sn) = 1
      code_new = code
      dnum    = 0
      do i = 1, norbs
          if(code_new(i) .eq. 1)then
              dnum = dnum + 2**(i-1)
          endif
      enddo

  end subroutine get_newstate

  subroutine get_power(n,bt,tn,m)
      implicit none
! which is the target number      
      integer,    intent(in)  ::  n

! which is the base
      integer,    intent(in)  ::  bn

! which is the top boundry of power
      integer,    intent(in)  ::  tn

! which is the output: poower
      integer,    intent(out) ::  m

      integer   :: i, tmp
      
      do i = 1, tn
          if(bn**i .eq. n)then
              m = i+1
          else
              stop "get power fail"
          endif
      enddo

  end subroutine get_power
!
  subroutine make_Cg(D1_dim,dp_num,dp_type,dp_vec)
      implicit none
! the dimension of every irrep of C4v
      integer,   intent(in)      :: D1_dim(5)

! this part will prepare the CG efficient of two irrep
! the number of reduction sub irreps when ith irrep and jth irrep direct product
      integer,   intent(out)     :: dp_num(5,5)

! the dimension of every reduction irrep , the number of reduction sub-irrep is 4 at most
      integer,   intent(out)     :: dp_type(5,5,4)

! the Cg Coefficients: the number of one direct product matrix; another direct
! product matrix; the number of reduction matrix; the row label; the column
! labels
      complex(8),intent(out)     :: dp_vec(5,5,4,4,2)

! constants
      u = 1.d0/sqrt(2.d0)

! local variables
      integer     :: i, j, k, ii, jj, kk, n, nn
      integer     :: basis1, basis2, serial1(norbs)
      complex     :: cone,czero,cima

! initial
      cone    = dcmplx(1.d0,0.d0)
      czero   = dcmplx(0.d0,0.d0)
      cima    = dcmplx(0.d0,1.d0)
      dp_num  = 0
      dp_type = 0
      dp_vec  = dcmplx(0.d0,0.d0)

! key part   
! A1-A1
      dp_num(1,1)         = 1
      dp_type(1,1,1)      = 1
      dp_vec(1,1,1,1,1)   = cone
! A1-A2
      dp_num(1,2)         = 1
      dp_type(1,2,1)      = 2
      dp_vec(1,2,1,1,1)   = cone
! A1-B1
      dp_num(1,3)         = 1
      dp_type(1,3,1)      = 3
      dp_vec(1,3,1,1,1)   = cone
! A1-B2
      dp_num(1,4)         = 1
      dp_type(1,4,1)      = 4
      dp_vec(1,4,1,1,1)   = cima
! A1-E
      dp_num(1,5)         = 1
      dp_type(1,5,1)      = 5
      dp_vec(1,5,1,1,1)   =-u
      dp_vec(1,5,1,1,2)   =-u
      dp_vec(1,5,1,2,1)   =-u*cima
      dp_vec(1,5,1,2,2)   = u*cima
! A2-A2
      dp_num(2,2)         = 1
      dp_type(2,2,1)      = 1
      dp_vec(2,2,1,1,1)   = cone
! A2-B1
      dp_num(2,3)         = 1
      dp_type(2,3,1)      = 4
      dp_vec(2,3,1,1,1)   = cone
! A2-B2
      dp_num(2,4)         = 1
      dp_type(2,4,1)      = 3
      dp_vec(2,4,1,1,1)   = cima
! A2-E
      dp_num(2,5)         = 1
      dp_type(2,5,1)      = 5
      dp_vec(2,5,1,1,1)   =-u
      dp_vec(2,5,1,1,2)   = u
      dp_vec(2,5,1,2,1)   =-u*cima
      dp_vec(2,5,1,2,2)   =-u*cima
! B1-B1
      dp_num(3,3)         = 1
      dp_type(3,3,1)      = 1
      dp_vec(3,3,1,1,1)   = cone
! B1-B2
      dp_num(3,4)         = 1
      dp_type(3,4,1)      = 2
      dp_vec(3,4,1,1,1)   = cima
! B1-E
      dp_num(3,5)         = 1
      dp_type(3,5,1)      = 5
      dp_vec(3,5,1,1,1)   =-u
      dp_vec(3,5,1,1,2)   =-u
      dp_vec(3,5,1,2,1)   = u*cima
      dp_vec(3,5,1,2,2)   =-u*cima
! B2-B2
      dp_num(4,4)         = 1
      dp_type(4,4,1)      = 1
      dp_vec(4,4,1,1,1)   =-cone
! B2-E
      dp_num(4,5)         = 1
      dp_type(4,5,1)      = 5
      dp_vec(4,5,1,1,1)   =-u*cima
      dp_vec(4,5,1,2,1)   =-u
      dp_vec(4,5,1,1,2)   = u*cima
      dp_vec(4,5,1,2,2)   =-u
! E-E
      dp_num(5,5)         = 4
      dp_type(5,5,1)      = 1
      dp_vec(5,5,1,1,1)   = u
      dp_vec(5,5,1,4,1)   = u
      dp_type(5,5,2)      = 2
      dp_vec(5,5,2,2,1)   =-u*cima
      dp_vec(5,5,2,3,1)   = u*cima
      dp_type(5,5,3)      = 3
      dp_vec(5,5,3,1,1)   = 0.5
      dp_vec(5,5,3,2,1)   = 0.5*cima
      dp_vec(5,5,3,3,1)   = 0.5*cima
      dp_vec(5,5,3,4,1)   =-0.5
      dp_type(5,5,4)      = 4
      dp_vec(5,5,4,1,1)   = 0.5
      dp_vec(5,5,4,2,1)   =-0.5*cima
      dp_vec(5,5,4,3,1)   =-0.5*cima
      dp_vec(5,5,4,4,1)   =-0.5
!-------------------------------------
! A2-A1
      dp_num(2,1)         = 1
      dp_type(2,1,1)      = 2
      dp_vec(2,1,1,1,1)   = cone
! B1-A1
      dp_num(3,1)         = 1
      dp_type(3,1,1)      = 3
      dp_vec(3,1,1,1,1)   = cone
! B2-A1
      dp_num(4,1)         = 1
      dp_type(4,1,1)      = 4
      dp_vec(4,1,1,1,1)   = cima
! E-A1
      dp_num(5,1)         = 1
      dp_type(5,1,1)      = 5
      dp_vec(5,1,1,1,1)   =-u
      dp_vec(5,1,1,1,2)   =-u
      dp_vec(5,1,1,2,1)   =-u*cima
      dp_vec(5,1,1,2,2)   = u*cima
! B1-A2
      dp_num(3,2)         = 1
      dp_type(3,2,1)      = 4
      dp_vec(3,2,1,1,1)   = cone
! B2-A2
      dp_num(4,2)         = 1
      dp_type(4,2,1)      = 3
      dp_vec(4,2,1,1,1)   = cima
! E-A2
      dp_num(5,2)         = 1
      dp_type(5,2,1)      = 5
      dp_vec(5,2,1,1,1)   =-u
      dp_vec(5,2,1,1,2)   = u
      dp_vec(5,2,1,2,1)   =-u*cima
      dp_vec(5,2,1,2,2)   =-u*cima
! B2-B1
      dp_num(4,3)         = 1
      dp_type(4,3,1)      = 2
      dp_vec(4,3,1,1,1)   = cima
! E-B1
      dp_num(5,3)         = 1
      dp_type(5,3,1)      = 5
      dp_vec(5,3,1,1,1)   =-u
      dp_vec(5,3,1,1,2)   =-u
      dp_vec(5,3,1,2,1)   = u*cima
      dp_vec(5,3,1,2,2)   =-u*cima
! E-B2
      dp_num(5,4)         = 1
      dp_type(5,4,1)      = 5
      dp_vec(5,4,1,1,1)   =-u*cima
      dp_vec(5,4,1,2,1)   =-u
      dp_vec(5,4,1,1,2)   = u*cima
      dp_vec(5,4,1,2,2)   =-u
  end subroutine make_Cg
