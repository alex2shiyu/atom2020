!>>>>author:sypeng@iphy.ac.cn
      subroutine gw_make_newUI2(num_orb,num_config,occ1,occ2,UNmtrxp00,basis,invcd,invsn,UImtrxp00)
      implicit none
      integer, external :: state_pick
      integer,      intent(in)   :: num_orb
      integer,      intent(in)   :: num_config
      integer,      intent(in)   :: occ1
      integer,      intent(in)   :: occ2
      integer,      intent(in)   :: basis(num_config)
      integer,      intent(in)   :: invcd(num_orb,num_config)
      integer,      intent(in)   :: invsn(0:2**num_orb-1)
      complex(8),   intent(in)   :: UNmtrxp00(num_orb,num_orb)
      complex(8),   intent(out)  :: UImtrxp00(num_config,num_config)
!
      integer:: i,j,n,m,num1,num,ii,jj,occt
      integer:: offset_n,offset_m,flag1,flag2,occ_n,cnt_flag,sgn1,n_fk,flag3(num_orb)
      integer:: code_a(Num_orb),code_old(Num_orb),code_new(Num_orb)
      integer:: code_last(6000,Num_orb),code_tmp(Num_orb)
      integer:: long1(Num_orb,Num_orb),lgth(num_orb)
      complex(8):: value_last(6000),v_tmp
!     
      UImtrxp00=dcmplx(0.d0,0.d0)
      lgth=0
      num1=0
      do n=occ1,occ2
         num1=num1+state_pick(n,num_orb)
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
         occt=occ_n
         num=state_pick(occ_n,num_orb)
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
	                     if(abs(UNmtrxp00(ii,j)).GT.0.00001)then
                             flag2=flag2+1
			                 long1(flag2,flag1)=ii
		                 end if
		            end do
		        else
		            cycle
                end if
		        lgth(flag1)=flag2
            end do!end every single wavefunction
            if(occt .eq. 7 .and. i .eq. num)then 
                write(20,*)'long1',long1
                do ii = 1, num_orb
                    do jj = 1, num_orb
                        write(22,'(3I4,2x)')jj,ii,long1(jj,ii)
                    enddo
                enddo
                write(15,*)'lgth',lgth
            endif
!           call code_tranp(occ,num_orb,long,lgth,code_last)
            call code_tranp(occt,occt,num_orb,long1,lgth,cnt_flag,code_last)
!	        print *,"-----------test5-below-code-tranp--------------"
            if(occt .eq. 7 .and. i .eq. num)then 
                write(31,*)'cnt_flag',cnt_flag
                do ii = 1, cnt_flag
                     write(32,'(14I2)')code_last(ii,:)
                enddo
            endif
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
		        call code_ord2(code_new,occt,Num_orb,num_config,basis,invcd,invsn,n_fk)
		        UImtrxp00(n_fk,i+offset_n)=v_tmp*sgn1+UImtrxp00(n_fk,i+offset_n)
	        end do
!!
!            open(23,file = 'code_syp.dat',access = 'append')
!            write(23,"(14I2)")code_a(:)
!            close(23)
            call code_shift3(code_a,num_config,occ_n,basis,invcd,invsn,Num_orb)
         end do!end n subspace
	     offset_n=offset_n+num
      end do !end n=occ1,occ2
      if(occ1.EQ.0)then
         UImtrxp00(1,1)=dcmplx(0.d0,0.d0)
      end if
      end subroutine gw_make_newUI2
!
      recursive subroutine code_tranp(occ,n,num_orb,long,lgth,cnt_flag,code_last)
!         
      implicit none
      integer,   intent(in)    ::  occ
      integer,   intent(in)    ::  n
      integer,   intent(in)    ::  num_orb
      integer,   intent(in)    ::  long(num_orb,num_orb)
      integer,   intent(in)    ::  lgth(num_orb)
      integer,   intent(out)   ::  cnt_flag
      integer,   intent(out)   ::  code_last(6000,num_orb)

      integer:: i,j,m,num,ii,jj
      integer:: offset_n,offset_m,flg1,flg2,occ_n
      integer:: code1_tmp(Num_orb),code_old(Num_orb)
      integer:: code_tmp(6000,Num_orb)

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
      call code_tranp(occ,n-1,num_orb,long,lgth,cnt_flag,code_last)
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

      subroutine code_shift3(code,num_config,occ_num,basis,invcd,invsn,Num_orb)
      implicit none
      integer  :: Num_orb,num_config
      integer  :: basis(num_config)
      integer  :: invcd(num_orb,num_config)
      integer  :: invsn(0:2**num_orb-1)
      integer:: occ_num,code(Num_orb),i,n,flag,flag1

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

      subroutine code_shift2(code,occ_num,Num_orb)
      implicit none
      integer:: Num_orb
      integer:: occ_num,code(Num_orb),i,n,flag,flag1

      flag=0
      flag1=1
      do i = 1, num_orb-1
          flag1 = 1
          if(code(i) .eq. 1)then
              do n=i,i+occ_num-1
                  flag1=flag1*code(n)
              enddo
              if(flag1 .eq. 1 .and. occ_num .gt. 1 .and. code(1) .eq. 0)then
                  code = 0
                  do n = 1,occ_num-1
                      code(n) = 1
                  enddo
                  code(i+occ_num) =1
                  flag=1
              else
                  do n = i+1,num_orb
                      if(code(n) .eq. 0)then
                          code(n) = 1
                          code(n-1) = 0
                          flag = 1
                          exit
                      endif
                  enddo
              endif
          endif
          if(flag .eq. 1)exit
      enddo

      end subroutine code_shift2

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCC this subroutine is very smart                                   CCC
      subroutine code_ord(code,occ_num,Num_orb,n)
      implicit none
      integer, external :: state_pick
      integer:: Num_orb,code(num_orb),occ_num,n,ii
      integer:: flg,nn,m,i,j,pos(0:occ_num)

      n=1
      i=0
      pos(0)=0
      do j=1,num_orb
          if(code(j).EQ.1)then
              i=i+1
	          pos(i)=j
	      end if
	  end do
	  do i=1,occ_num
          if(pos(i).GT.(pos(i-1)+1))then
              do ii=pos(i-1)+1,pos(i)-1
                  n=n+state_pick(occ_num-i,num_orb-ii)
	          end do
	      end if
	  end do

      end subroutine code_ord


      subroutine code_ord2(code,occ_num,Num_orb,num_config,basis,invcd,invsn,n)
      implicit none
      integer, external :: state_pick
      integer  :: Num_orb,num_config
      integer  :: basis(num_config)
      integer  :: invcd(num_orb,num_config)
      integer  :: invsn(0:2**num_orb-1)
      integer  :: code(num_orb),occ_num,n,ii
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

      subroutine code_shift(code,occ_num,Num_orb)
      implicit none
      integer:: Num_orb
      integer:: occ_num,code(Num_orb),i,n,flag,index,flag1


      index=0
      flag=0

      do n=Num_orb,1,-1
         if(code(n).eq.1)then
            index=index+1
            if((n+index.le.Num_orb).and.(flag.eq.0))then
                 code(n)=0
                 code(n+1)=1
                 if(index.gt.1)then
                 do i=n+2,n+index
                    code(i)=1
                 end do
                 do i=n+index+1,Num_orb
                    code(i)=0
                 end do
                 end if 

                             
                
                 flag=1
            end if
          end if

       end do

       end subroutine code_shift

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
