subroutine trans_code(code,Num_orb,n)
    implicit none
    integer, intent(in)  :: Num_orb
    integer, intent(in)  :: code(Num_orb)
    integer, intent(out) :: n
    integer :: i, occ_num

!f2py intent(in)  num_orb
!f2py intent(in)  code
!f2py intent(out) n
!f2py depend(num_orb) code

    n=0
    do i=0,num_orb-1
        if(code(i+1).eq.1)n=n+2**i
    end do
end subroutine trans_code

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine inverse_trans_code(code,Num_orb,nn)
    implicit none
    integer,  intent(in)  :: Num_orb
    integer,  intent(out) :: code(Num_orb)
    integer,  intent(in)  :: nn
    integer   ::   i, n, n1

!f2py  intent(in)   num_orb
!f2py  intent(in)   nn
!f2py  intent(out)  code
!f2py  depend(num_orb) code


	n=nn

    do i=1,num_orb
        n1=n/2
        code(i)=n-n1*2
        n=n1
    end do
end subroutine inverse_trans_code

 

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine code_shift(code,occ_num,Num_orb)
    implicit none
    integer,  intent(in)     :: Num_orb
    integer,  intent(in)     :: occ_num
    integer,  intent(inout)  :: code(Num_orb)
    integer   ::  i,n,flag,index1

!f2py  intent(in)   num_orb    
!f2py  intent(in)   occ_num    
!f2py  intent(inout)   code    
!f2py  depend(num_orb) code
    index1=0
    flag=0

    do n=Num_orb,1,-1
       if(code(n).eq.1)then
       index1=index1+1
          if((n+index1.le.Num_orb).and.(flag.eq.0))then
               code(n)=0
               code(n+1)=1
               if(index1.gt.1)then
               do i=n+2,n+index1
                  code(i)=1
               end do
               do i=n+index1+1,Num_orb
                  code(i)=0
               end do
               end if 
               flag=1
          end if
        end if
     end do
 end subroutine code_shift

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCC                                                                 CCC



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       
subroutine code_comp(code_alpha,code_beta,Num_orb,res)
    implicit none
 
    integer,  intent(in)  :: Num_orb
    integer,  intent(in)  :: code_alpha(Num_orb)
    integer,  intent(in)  :: code_beta(Num_orb)
    integer,  intent(out) :: res
    integer:: k

!f2py  intent(in)    num_orb
!f2py  intent(in)    code_alpha
!f2py  intent(in)    code_beta
!f2py  intent(out)   res
!f2py  depend(num_orb) code_alpha
!f2py  depend(num_orb) code_beta

    res=1  

    do k=1,Num_orb
       if(code_beta(k).ne.code_alpha(k))res=0
    end do

end subroutine code_comp
            

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCC                                                                 CCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       
subroutine code_action2(code_alpha,code_beta,Num_orb,n1,n2,m1,m2,res)
    implicit none
    integer,   intent(in)  :: num_orb
    integer,   intent(in)  :: code_alpha(Num_orb)
    integer,   intent(in)  :: n1, n2, m1, m2
    integer,   intent(out) :: res
    integer,   intent(out) :: code_beta(Num_orb)
    integer   :: k,sgn,i,r1,r2

!f2py intent(in)    num_orb    
!f2py intent(in)    n1    
!f2py intent(in)    n2    
!f2py intent(in)    m1   
!f2py intent(in)    m2  
!f2py intent(in)    code_alpha    
!f2py intent(out)   res    
!f2py intent(out)   code_beta    
!f2py depend(num_orb)  code_alpha
!f2py depend(num_orb)  code_beta

    sgn=1
    res=0
    r1=0
    r2=0
    
    code_beta=code_alpha

    if((code_beta(m1).eq.1).and.(code_beta(m2).eq.1))then
       r1=1
       do i=1,m1-1
          if(code_beta(i).eq.1)r1=r1
       end do
       code_beta(m1)=0

       do i=1,m2-1
          if(code_beta(i).eq.1)r1=r1
       end do
       code_beta(m2)=0


    if((code_beta(n1).eq.0).and.(code_beta(n2).eq.0))then

       r2=1
       do i=1,n2-1
          if(code_beta(i).eq.1)r2=r2
       end do
       code_beta(n2)=1

       do i=1,n1-1
          if(code_beta(i).eq.1)r2=r2
       end do
       code_beta(n1)=1

     end if
    
     end if

     res=r1*r2

end subroutine code_action2


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCC                                                                 CCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCC                                                                 CCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       
subroutine code_action1(code_alpha,code_beta,Num_orb,n1,res)
    implicit none
    integer,   intent(in)  :: num_orb
    integer,   intent(in)  :: code_alpha(Num_orb)
    integer,   intent(out) :: code_beta(Num_orb)
    integer,   intent(in)  :: n1
    integer,   intent(out) :: res
    integer    :: r2
!f2py intent(in)   num_orb
!f2py intent(in)   code_alpha 
!f2py intent(in)   n1
!f2py intent(out)  code_beta
!f2py intent(out)  res
!f2py depend(num_orb) code_alpha
!f2py depend(num_orb) code_beta


    res=0
    r2=0
    
    code_beta=code_alpha
    if(code_beta(n1).eq.0)then
       r2=1
       code_beta(n1)=1
    end if
    

    res=r2

end subroutine code_action1


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCC                                                                 CCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine sym_array(num_orb,array,sym)
    implicit none

	integer,  intent(in)    :: num_orb
    integer,  intent(in)    :: sym(num_orb)
    real*8,   intent(inout) :: array(num_orb)
	integer:: i,j,flag,p1,p2
	real*8::  new,r,res

!f2py intent(in)       num_orb
!f2py intent(in)       sym
!f2py intent(inout)    array
!f2py depend(num_orb)  sym
!f2py depend(num_orb)  array

    do i=1,num_orb
        if(abs(array(i)).lt.0.01)array(i)=0.01
	end do
    p1=1
111	flag=sym(p1)
    res=0
	p2=p1
	
	do i=p1,num_orb
	   if(sym(i).eq.flag)then
           res=res+array(i)
           p2=i
       endif
	end do

	new=res/real(p2-p1+1)
	do i=p1,p2
	   array(i)=new
	end do

	p1=p2+1
	if(p1.le.num_orb)goto 111

end subroutine sym_array

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine M0func(norb,ncfg,Cg,n0vec,m0fock)
    implicit none
    integer,  intent(in)  :: norb
    integer,  intent(in)  :: ncfg  
    integer,  intent(in)  :: Cg(ncfg)
    real(kind=8),  intent(in) , dimension(norb)  :: n0vec
    real(kind=8),  intent(out) , dimension(ncfg) :: m0fock

!f2py intent(in)    norb
!f2py intent(in)    ncfg
!f2py intent(in)      Cg
!f2py intent(in)   n0vec
!f2py intent(out) m0fock  
!f2py depend(ncfg)    cg
!f2py depend(norb) n0vec
!f2py depend(ncfg) m0fock

! local variables:
    integer :: n,m,i,j,alpha,beta,i1,fsgn,nn
    integer :: statealpha(norb)
    real(kind=8)  :: res
    do n=1,ncfg
        nn=Cg(n)
        call inverse_trans_code(statealpha,norb,nn)
        res=1.D0
        do i=1,norb
            if(statealpha(i).eq.0)then
                res=res*(1.D0-n0vec(i))
            else
                res=res*n0vec(i)
            endif
        enddo
        m0fock(n)=res
    enddo
end subroutine M0func
 
subroutine dMdN0func(norb,ncfg,Cg,n0vec,m0fock,dMdN0)
    implicit none
    integer,  intent(in)  :: norb
    integer,  intent(in)  :: ncfg  
    integer,  intent(in)  :: Cg(ncfg)
    real(kind=8),  intent(in) , dimension(norb)  :: n0vec
    real(kind=8),  intent(in) , dimension(ncfg) :: m0fock
    real(kind=8),  intent(out)  :: dMdN0(ncfg,norb)

!f2py intent(in)     norb
!f2py intent(in)     ncfg
!f2py intent(in)       Cg
!f2py intent(in)    n0vec
!f2py intent(in)   m0fock  
!f2py intent(out)  dMdN0V  
!f2py depend(ncfg)     cg
!f2py depend(norb)  n0vec
!f2py depend(ncfg) m0fock
!f2py depend(ncfg,norb) dMdN0V

! local variables:
    integer :: n,m,i,j,alpha,beta,i1,fsgn,nn
    integer :: statealpha(norb)
    real(kind=8)  :: res
    do n=1,ncfg
        nn=Cg(n)
        call inverse_trans_code(statealpha,norb,nn)
        res=1.D0
        do i=1,norb
            if(statealpha(i).eq.0)then
                res=m0fock(n)/(n0vec(i) - 1.0)
            else
                res=m0fock(n)/n0vec(i)
            endif
            dMdN0(n,i)=res
        enddo
    enddo
end subroutine dMdN0func
 


subroutine Soptfunc(num_orb,num_config,Cg,fp_add,fp_remove)
    implicit none
    integer , intent(in)  :: num_orb
    integer , intent(in)  :: num_config  
    integer , intent(in)  :: Cg(num_config)  
    integer , intent(out) :: fp_add(num_orb,num_config)
    integer , intent(out) :: fp_remove(num_orb,num_config)
    integer :: alpha,beta,nn,n,m,i,j,fsgn
    integer :: state_alpha(num_orb),state_beta(num_orb)

!f2py intent(in)  num_orb
!f2py intent(in)  num_config
!f2py intent(in)  Cg
!f2py intent(out) fp_add
!f2py intent(out) fp_remove
!f2py depend(num_config) Cg
!f2py depend(num_orb,num_config) fp_add
!f2py depend(num_orb,num_config) fp_remove
    fp_add=0
    fp_remove=0
    do alpha=1,num_config
        nn=Cg(alpha)
    	call inverse_trans_code(state_alpha,num_orb,nn)
        do n=1,num_orb
            if(state_alpha(n).eq.0)then
                state_beta=state_alpha
		        fsgn=1
		        do m=1,n-1
                    if(state_beta(m).eq.1)fsgn=-fsgn
		        enddo
		        state_beta(n)=1
		        call trans_code(state_beta,Num_orb,nn)
		        do m=1,num_config
                    if(Cg(m).eq.nn)then
                        fp_add(n,alpha)=m*fsgn
                        exit
		            endif
		        enddo
            else
		        state_beta=state_alpha
		        fsgn=1
		        do m=1,n-1
                    if(state_beta(m).eq.1)fsgn=-fsgn
		        enddo
		        state_beta(n)=0
		        call trans_code(state_beta,Num_orb,nn)
		        do m=1,num_config
                    if(Cg(m).eq.nn)then
                        fp_remove(n,alpha)=m*fsgn
                        exit
                    endif
                enddo
            endif
        enddo
    enddo
    return
end subroutine Soptfunc
