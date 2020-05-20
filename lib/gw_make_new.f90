!>>>>author:sypeng@iphy.ac.cn
      subroutine gw_make_newui(num_orb,num_config,totncfgs,unitary_len,occ1,occ2,UNmtrxp00,basis,invcd,invsn,prec,UImtrxp00)
      implicit none
!     integer, external :: state_pick
      integer,      intent(in)   :: num_orb
      integer,      intent(in)   :: num_config
      integer,      intent(in)   :: totncfgs
      integer,      intent(in)   :: unitary_len
      integer,      intent(in)   :: occ1
      integer,      intent(in)   :: occ2
      integer,      intent(in)   :: basis(num_config)
      integer,      intent(in)   :: invcd(num_orb,num_config)
      integer,      intent(in)   :: invsn(0:totncfgs-1)
      complex(8),   intent(in)   :: UNmtrxp00(num_orb,num_orb)
      real(kind=8), optional     :: prec
      complex(8),   intent(out)  :: UImtrxp00(num_config,num_config)
!
      integer:: i,j,n,m,num1,num,ii,jj,occt,int_tmp
      integer:: offset_n,offset_m,flag1,flag2,occ_n,cnt_flag,sgn1,n_fk,flag3(num_orb)
      integer:: code_a(Num_orb),code_old(Num_orb),code_new(Num_orb)
      integer:: code_last(unitary_len,Num_orb),code_tmp(Num_orb)
      integer:: long1(Num_orb,Num_orb),lgth(num_orb)
      complex(8):: value_last(unitary_len),v_tmp
      real(kind=8) :: prec_

!f2py intent(in) num_orb
!f2py intent(in) num_config
!f2py intent(in) totncfgs
!f2py intent(in) unitary_len
!f2py intent(in) occ1
!f2py intent(in) occ2
!f2py intent(in) basis
!f2py intent(in) invcd
!f2py intent(in) invsn
!f2py intent(in) prec
!f2py intent(in) UNmtrxp00
!f2py intent(out) UImtrxp00
!f2py depend(num_config) basis
!f2py depend(num_orb,num_config) invcd
!f2py depend(totncfgs) invsn
!f2py depend(num_orb)  UNmtrxp00
!f2py depend(num_config)  UImtrxp00
      !     
      UImtrxp00=dcmplx(0.d0,0.d0)
      lgth=0
      num1=0
      prec_ = 1.0E-5
      if(present(prec)) prec_ = prec
      do n=occ1,occ2
         call state_pick(n,num_orb,int_tmp)
!        num1=num1+state_pick(n,num_orb)
         num1=num1+int_tmp
      end do
!
      if(num_config.NE.num1)then
          print *,"Running error in gw_make_new(:newUI2),error code 001"
          stop
      end if
      offset_n=0
      offset_m=0
!      
      do occ_n=occ1,occ2
         print *,'======================================'
         occt=occ_n
         call state_pick(occ_n,num_orb,int_tmp)
!        num=state_pick(occ_n,num_orb)
         num=int_tmp
    	 code_last=0
         do i=1,num_orb
             if(i.le.occ_n)then
                 code_a(i)=1
	         else
	             code_a(i)=0
	         end if
	     end do
!         
         do i=1,num
	        long1=0
	        lgth=0
	        flag1=0
	        flag3=0
	        cnt_flag=0
            do j=1,num_orb
                if(code_a(j).EQ.1)then
                    flag1=flag1+1
		            flag3(flag1)=j
		            flag2=0
                    do ii=1,num_orb
!              UNmtrx must be obtained by eigen function
	                     if(abs(UNmtrxp00(ii,j)).GT. prec_ )then
                             flag2=flag2+1
			                 long1(flag2,flag1)=ii
		                 end if
		            end do
		        else
		            cycle
                end if
		        lgth(flag1)=flag2
            end do!end every single wavefunction
!!          if(occt .eq. 7 .and. i .eq. num)then 
!!              write(20,*)'long1',long1
!!              do ii = 1, num_orb
!!                  do jj = 1, num_orb
!!                      write(22,'(3I4,2x)')jj,ii,long1(jj,ii)
!!                  enddo
!!              enddo
!!              write(15,*)'lgth',lgth
!!          endif
!           call code_tranp(occ,num_orb,long,lgth,code_last)
            print *, '-------------------'
            print *,'nmin-nmax:',occ_n,'ncfgs:',i
            call code_tranp(occt,occt,num_orb,unitary_len,long1,lgth,cnt_flag,code_last)
            print *, ''
            print *, ''
!	        print *,"-----------test5-below-code-tranp--------------"
!!          if(occt .eq. 7 .and. i .eq. num)then 
!!              write(31,*)'cnt_flag',cnt_flag
!!              do ii = 1, cnt_flag
!!                   write(32,'(14I2)')code_last(ii,:)
!!              enddo
!!          endif
!!
            do ii=1,cnt_flag
                code_tmp=0
		        code_new=0
		        v_tmp=dcmplx(1.d0,0.d0)
!
	            do j=1,occt
                    code_tmp(j)=code_last(ii,j)
		            v_tmp=v_tmp*UNmtrxp00(code_last(ii,j),flag3(j))
		        end do
		        call sgnp(occt,num_orb,code_tmp,sgn1)
!
		        do j=1,occt
		            code_new(code_last(ii,j))=1
		        end do
!
		        call code_ord2(occt,Num_orb,num_config,totncfgs,code_new,basis,invcd,invsn,n_fk)
		        UImtrxp00(n_fk,i+offset_n)=v_tmp*sgn1+UImtrxp00(n_fk,i+offset_n)
	        end do
!!
!            open(23,file = 'code_syp.dat',access = 'append')
!            write(23,"(14I2)")code_a(:)
!            close(23)
            call code_shift3(num_config,totncfgs,occ_n,basis,invcd,invsn,Num_orb,code_a)
         end do!end n subspace
	     offset_n=offset_n+num
      end do !end n=occ1,occ2
      if(occ1.EQ.0)then
         UImtrxp00(1,1)=dcmplx(0.d0,0.d0)
      end if
      end subroutine gw_make_newui
!
      recursive subroutine code_tranp(occ,n,num_orb,unitary_len,long,lgth,cnt_flag,code_last)
!         
      implicit none
      integer,   intent(in)    ::  occ
      integer,   intent(in)    ::  n
      integer,   intent(in)    ::  num_orb
      integer,   intent(in)    ::  unitary_len
      integer,   intent(in)    ::  long(num_orb,num_orb)
      integer,   intent(in)    ::  lgth(num_orb)
      integer,   intent(out)   ::  cnt_flag
      integer,   intent(out)   ::  code_last(unitary_len,num_orb)

      integer:: i,j,m,num,ii,jj
      integer:: offset_n,offset_m,flg1,flg2,occ_n
      integer:: code1_tmp(Num_orb),code_old(Num_orb)
      integer:: code_tmp(unitary_len,Num_orb)

!      real*8:: U_mtrx(num_orb,num_orb,num_orb,num_orb)
!      complex(dp):: Elm_mtrx1(num_orb,num_orb)
!   
!      cnt_flag=0
!      cnt_flag=0
      flg2=0
      code_tmp=0
      if(n.EQ.0)then
          flg2=0
	      cnt_flag=flg2+cnt_flag
	      return
      else if(n.LT.0)then
          print *,"there is something wrong in code_tranp(1)!"
      else if(n.GT.occ)then
          print *,"there is something wrong in code_tranp(2)!"
      else if(n.EQ.occ)then
          flg2=0
          do i=1,lgth(occ)
             code_tmp(i,occ)=long(i,occ)
	         flg2=flg2+1
          end do
	      code_last=code_tmp
	      cnt_flag=flg2
      else 
          flg1=1
          flg2=0
          m=occ-n
          do ii=1,lgth(n)
              do jj=1,cnt_flag
                  flg1=1
                  do i=n+1,occ
                      if(long(ii,n).EQ.code_last(jj,i))then
                          flg1=0
           	          else
           	              continue
           	          end if
                  end do
                  if(flg1.EQ.0)then
                      cycle
                  else
                      flg2=flg2+1
                      do j=n+1,occ
                          code_tmp(flg2,j)=code_last(jj,j)
                      end do
                      code_tmp(flg2,n)=long(ii,n)
                  end if
              end do
          end do
	      code_last=code_tmp
	      cnt_flag=flg2
      end if 
      print *, 'ORBth = ',n,'cnt_flag=',cnt_flag
      call code_tranp(occ,n-1,num_orb,unitary_len,long,lgth,cnt_flag,code_last)
      end subroutine code_tranp



      subroutine sgnp(occ,num_orb,code1,res)
      implicit none
      integer, intent(in)   :: occ
      integer, intent(in)   :: num_orb
      integer, intent(in)   :: code1(num_orb)
      integer, intent(out)  :: res

      integer:: i,j,n,m,num,ii,jj,flag
      integer:: code_old(Num_orb),lgth(num_orb)
      
      flag=1
      if(occ.EQ.1)then
         flag=1
      else if(occ.EQ.0)then
         flag=1
      else
         do i=1,occ-1
             do j=i+1,occ
                 if(code1(i).GT.code1(j))then
                     flag=flag*(-1)
                 endif
             end do
         end do
      end if
!
!
      do i=1,occ
         if(code1(i).EQ.0)then
             print *,"there is mistake(1) in sgnp"
         end if
      end do
      do i=occ+1,num_orb
          if(code1(i).NE.0)then
              print *,"there is mistake(2) in sgnp"
          end if
      end do

      res=flag
      end subroutine sgnp

      subroutine code_shift3(num_config,totncfgs,occ_num,basis,invcd,invsn,Num_orb,code)
      implicit none
      integer, intent(in)    :: Num_orb,num_config,totncfgs
      integer, intent(in)    :: basis(num_config)
      integer, intent(in)    :: invcd(num_orb,num_config)
      integer, intent(in)    :: invsn(0:totncfgs-1)
      integer, intent(in)    :: occ_num
      integer, intent(inout) :: code(Num_orb)
      
      integer :: i,n,flag,flag1

      flag = 0
      do i = 1,num_orb
          if(code(i) .eq. 1)then
              flag = flag+2**(i-1)
          endif
      enddo
      do i = flag+1 , 2**num_orb-1
          call verif( i, num_orb, n )
          if (n .eq. occ_num)then
              code = invcd(:,invsn(i))
              exit
          endif
      enddo
      return
      end subroutine code_shift3


      subroutine code_ord2(occ_num,Num_orb,num_config,totncfgs,code,basis,invcd,invsn,n)
      implicit none
!     integer, external :: state_pick
      integer, intent(in)  :: Num_orb
      integer, intent(in)  :: num_config
      integer, intent(in)  :: totncfgs
      integer, intent(in)  :: code(num_orb)
      integer, intent(in)  :: occ_num
      integer, intent(in)  :: basis(num_config)
      integer, intent(in)  :: invcd(num_orb,num_config)
      integer, intent(in)  :: invsn(0:totncfgs-1)
      integer, intent(out) :: n 
      !
      integer  :: ii
      integer  :: flg,nn,m,i,j,pos(0:occ_num)

      nn=0
      m=0
      do ii = 1, num_orb
          if(code(ii) .eq. 1)then
              nn = nn + 2**(ii-1)
          endif
      enddo
      n=invsn(nn)

      end subroutine code_ord2

!>>> calculate combination algebra 
  subroutine state_pick(ntiny, nlarg,value1)
     implicit none

! external variables
     integer, intent(in)  :: ntiny
     integer, intent(in)  :: nlarg
     integer, intent(out) :: value1

! local variables
     integer :: i

! auxiliary integer variable
     integer :: nlow

! numberator of the combination algebra
     real(8) :: numer

! denominator of the combination algebra
     real(8) :: denom

! result value of the combination algebra
     integer :: value

! transform the combination algebra
     nlow = min(ntiny, nlarg-ntiny)

! numerator in combination algebra
     numer = 1.d0
     do i=nlarg-nlow+1,nlarg
        numer = numer * dble(i)
     enddo ! over i={nlarg-nlow+1,nlarg} loop

! denominator in combination algebra
     denom = 1.d0
     do i=1,nlow
        denom = denom * dble(i)
     enddo ! over i={1,nlow} loop

! result value
     value1 = nint(numer / denom)

     return
  end subroutine  state_pick

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

! number of total electrons
     nbits = 0
     do ipos=1,norbs
         if( btest(inum, ipos-1) ) nbits = nbits + 1
     enddo ! over ipos={1,norbs} loop

     return
  end subroutine verif
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCC this subroutine is very smart                                   CCC
!      subroutine code_ord(code,occ_num,Num_orb,n)
!      implicit none
!      integer, external :: state_pick
!      integer:: Num_orb,code(num_orb),occ_num,n,ii
!      integer:: flg,nn,m,i,j,pos(0:occ_num)
!
!      n=1
!      i=0
!      pos(0)=0
!      do j=1,num_orb
!          if(code(j).EQ.1)then
!              i=i+1
!	          pos(i)=j
!	      end if
!	  end do
!	  do i=1,occ_num
!          if(pos(i).GT.(pos(i-1)+1))then
!              do ii=pos(i-1)+1,pos(i)-1
!                  n=n+state_pick(occ_num-i,num_orb-ii)
!	          end do
!	      end if
!	  end do
!
!      end subroutine code_ord
