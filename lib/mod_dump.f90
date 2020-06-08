    subroutine dump_1d2dc(a,b,mat,path,nline,ioflag,prec)
        implicit none
! first nline lines which are introduction not data
        integer,    optional    :: nline

! the dimension of 3d matrix
        integer,    intent(in)  :: a
                 
        integer,    intent(in)  :: b


! the dump matrix
        complex(kind=8), intent(in) :: mat(a*b)

! the dump path
        character(len = *),  intent(in) :: path

! io flag: 1:replace, 2:append
        integer,    optional    :: ioflag

! control prec of output
        real(kind=8),   optional    :: prec

! local variables
        integer :: iatm, ibase1, ibase2, iorb, jorb, korb
        integer :: ioflag_, nline_
        real(kind=8) :: prec_
!f2py intent(in) a
!f2py intent(in) b
!f2py intent(in) nline
!f2py intent(in) path
!f2py intent(in) ioflag
!f2py intent(in) prec
!f2py intent(in) mat
!f2py depend(a,b) mat

! initial
        ioflag_ = 1
        nline_  = 0
        prec_   = 1.0d-16
        if(present(ioflag)) ioflag_ = ioflag
        if(present(nline))  nline_  = nline
        if(present(prec))   prec_   = prec
    
        
! key part
        if(ioflag_ .eq. 2)then
            open(40,file = path, access="append")
        elseif(ioflag_ .eq. 1)then
            open(40,file = path, status="replace")
        else
            open(40,file = path, status="replace")
        endif
        if(nline_ .gt. 0)then
            do iorb =1, nline_
                write(40,*) ''
            enddo
        else
            continue
        endif
        do iatm= 1,a
            ibase1 = (iatm-1)*b
            do iorb = 1,b
                if(abs(mat(ibase1+iorb)) .gt. prec_)then
                    write(40,'(2i10,2(7X,f35.16))')iatm,iorb,real(mat(ibase1+iorb)),aimag(mat(ibase1+iorb))
                endif
            enddo
        enddo
        close(40)

    end subroutine dump_1d2dc
    
    subroutine dump_1d2dr(a,b,mat,path,nline,ioflag,prec)
        implicit none
! first nline lines which are introduction not data
        integer,    optional    :: nline

! the dimension of 3d matrix
        integer,    intent(in)  :: a
                 
        integer,    intent(in)  :: b

! the dump matrix
        real(kind=8), intent(in) :: mat(a*b)

! the dump path
        character(len = *),  intent(in) :: path

! io flag: 1:replace, 2:append
        integer,    optional      :: ioflag

! control prec of output
        real(kind=8),   optional    :: prec

! local variables
        integer :: iatm, ibase1, ibase2, iorb, jorb, korb
        integer :: ioflag_, nline_
        real(kind=8) :: prec_
!f2py intent(in) a
!f2py intent(in) b
!f2py intent(in) nline
!f2py intent(in) path
!f2py intent(in) ioflag
!f2py intent(in) prec
!f2py intent(in) mat
!f2py depend(a,b) mat

        
! initial
        ioflag_ = 1
        nline_  = 0
        prec_   = 1.0d-16
        if(present(ioflag)) ioflag_ = ioflag
        if(present(nline))  nline_  = nline
        if(present(prec))   prec_   = prec
        
        if(ioflag_ .eq. 2)then
            open(40,file = path, access="append")
        elseif(ioflag_ .eq. 1)then
            open(40,file = path, status="replace")
        else
            open(40,file = path, status="replace")
        endif
        if(nline_ .gt. 0)then
            do iorb =1, nline_
                write(40,*) ''
            enddo
        else
            continue
        endif
        do iatm= 1,a
            ibase1 = (iatm-1)*b
            do iorb = 1,b
                if(abs(mat(ibase1+iorb)) .gt. prec_)then
                    write(40,'(2i10,7X,f35.16)')iatm,iorb,mat(ibase1+iorb)
                endif
            enddo
        enddo
        close(40)

    end subroutine dump_1d2dr
    !
    !
    subroutine dump_2d3dc(a,b,c,mat,path,nline,ioflag,prec)
        implicit none
! first nline lines which are introduction not data
        integer,    optional    :: nline

! the dimension of 3d matrix
        integer,    intent(in)  :: a
                 
        integer,    intent(in)  :: b
                 
        integer,    intent(in)  :: c


! the dump matrix
        complex(kind=8), intent(in) :: mat(a*b,a*c)

! the dump path
        character(len = *),  intent(in) :: path

! io flag: 1:replace, 2:append
        integer,    optional    :: ioflag

! control prec of output
        real(kind=8),   optional    :: prec

! local variables
        integer :: iatm, ibase1, ibase2, iorb, jorb, korb
        integer :: ioflag_, nline_
        real(kind=8) :: prec_
!f2py intent(in) a
!f2py intent(in) b
!f2py intent(in) c
!f2py intent(in) nline
!f2py intent(in) path
!f2py intent(in) ioflag
!f2py intent(in) prec
!f2py intent(in) mat
!f2py depend(a,b,c) mat


! initial
        ioflag_ = 1
        nline_  = 0
        prec_   = 1.0d-16
        if(present(ioflag)) ioflag_ = ioflag
        if(present(nline))  nline_  = nline
        if(present(prec))   prec_   = prec
        
! key part
        if(ioflag_ .eq. 2)then
            open(40,file = path, access="append")
        elseif(ioflag_ .eq. 1)then
            open(40,file = path, status="replace")
        else
            open(40,file = path, status="replace")
        endif
        if(nline_ .gt. 0)then
            do iorb =1, nline_
                write(40,*) ''
            enddo
        else
            continue
        endif
        do iatm= 1,a
            ibase1 = (iatm-1)*b
            ibase2 = (iatm-1)*c
            do iorb = 1,b
                do jorb = 1,c
                    if(abs(mat(ibase1+iorb,ibase2+jorb)) .gt. prec_)then
                        write(40,'(3i10,2(7X,f35.16))')iatm,iorb,jorb,&
                            real(mat(ibase1+iorb,ibase2+jorb)),&
                            aimag(mat(ibase1+iorb,ibase2+jorb))
                    endif
                enddo
            enddo
        enddo
        close(40)

    end subroutine dump_2d3dc
    
    subroutine dump_2d3dr(a,b,c,mat,path,nline,ioflag,prec)
        implicit none
! first nline lines which are introduction not data
        integer,    optional    :: nline

! the dimension of 3d matrix
        integer,    intent(in)  :: a
                 
        integer,    intent(in)  :: b
                 
        integer,    intent(in)  :: c


! the dump matrix
        real(kind=8), intent(in) :: mat(a*b,a*c)

! the dump path
        character(len = *),  intent(in) :: path

! io flag: 1:replace, 2:append
        integer,    optional      :: ioflag

! control prec of output
        real(kind=8),   optional    :: prec


! local variables
        integer :: iatm, ibase1, ibase2, iorb, jorb, korb
        integer :: ioflag_, nline_
        real(kind=8) :: prec_
!f2py intent(in) a
!f2py intent(in) b
!f2py intent(in) c
!f2py intent(in) nline
!f2py intent(in) path
!f2py intent(in) ioflag
!f2py intent(in) prec
!f2py intent(in) mat
!f2py depend(a,b,c) mat


        
! initial
        ioflag_ = 1
        nline_  = 0
        prec_   = 1.0d-16
        if(present(ioflag)) ioflag_ = ioflag
        if(present(nline))  nline_  = nline
        if(present(prec))   prec_   = prec
        
        if(ioflag_ .eq. 2)then
            open(40,file = path, access="append")
        elseif(ioflag_ .eq. 1)then
            open(40,file = path, status="replace")
        else
            open(40,file = path, status="replace")
        endif
        if(nline_ .gt. 0)then
            do iorb =1, nline_
                write(40,*) ''
            enddo
        else
            continue
        endif
        do iatm= 1,a
            ibase1 = (iatm-1)*b
            ibase2 = (iatm-1)*c
            do iorb = 1,b
                do jorb = 1,c
                    if(abs(mat(ibase1+iorb,ibase2+jorb)) .gt. prec_)then
                        write(40,'(3i10,7X,f35.16)')iatm,iorb,jorb,mat(ibase1+iorb,ibase2+jorb)
                    endif
                enddo
            enddo
        enddo
        close(40)

    end subroutine dump_2d3dr

  subroutine dump_4dc(a,b,c,d,mat,path,nline,ioflag,prec)
      implicit none
! the dimension of 3d matrix
     integer,    intent(in)  :: a
                 
     integer,    intent(in)  :: b
                 
     integer,    intent(in)  :: c

     integer,    intent(in)  :: d

! the dump matrix
     complex(kind=8), intent(in) :: mat(a,b,c,d)

! the dump path
     character(len = *),  intent(in) :: path

! first nline lines which are introduction not data
        integer,    optional    :: nline

! io flag: 1:replace, 2:append
        integer,    optional    :: ioflag

! control prec of output
        real(kind=8),   optional    :: prec


! local variables
        integer :: iatm, ibase1, ibase2, iorb, jorb, korb
        integer :: ioflag_, nline_
        real(kind=8) :: prec_
!f2py intent(in) a
!f2py intent(in) b
!f2py intent(in) c
!f2py intent(in) d
!f2py intent(in) nline
!f2py intent(in) path
!f2py intent(in) ioflag
!f2py intent(in) prec
!f2py intent(in) mat
!f2py depend(a,b,c,d) mat



! initial
        ioflag_ = 1
        nline_  = 0
        prec_   = 1.0d-16
        if(present(ioflag)) ioflag_ = ioflag
        if(present(nline))  nline_  = nline
        if(present(prec))   prec_   = prec
        
! key part
        if(ioflag_ .eq. 2)then
            open(40,file = path, access="append")
        elseif(ioflag_ .eq. 1)then
            open(40,file = path, status="replace")
        else
            open(40,file = path, status="replace")
        endif
        if(nline_ .gt. 0)then
            do iorb =1, nline_
                write(40,*) ''
            enddo
        else
            continue
        endif

       do iatm= 1,a
           do iorb = 1,b
               do jorb = 1,c
                   do korb = 1,d
                       if(abs(mat(iatm,iorb,jorb,korb)) .gt. prec_)then
                           write(40,'(4i10,2(7X,f35.16))')iatm,iorb,jorb,&
                               korb,real(mat(iatm,iorb,jorb,korb)),&
                               aimag(mat(iatm,iorb,jorb,korb))
                       endif
                   enddo
               enddo
           enddo
       enddo
     close(40)


  end subroutine dump_4dc

  subroutine dump_5dc(a,b,c,d,e,mat,path,nline,ioflag,prec)
      implicit none
! the dimension of 3d matrix
     integer,    intent(in)  :: a
                 
     integer,    intent(in)  :: b
                 
     integer,    intent(in)  :: c

     integer,    intent(in)  :: d

     integer,    intent(in)  :: e

! the dump matrix
     complex(kind=8), intent(in) :: mat(a,b,c,d,e)

! the dump path
     character(len = *),  intent(in) :: path

! first nline lines which are introduction not data
        integer,    optional    :: nline

! io flag: 1:replace, 2:append
        integer,    optional    :: ioflag

! control prec of output
        real(kind=8),   optional    :: prec


! local variables
        integer :: iatm, ibase1, ibase2, iorb, jorb, korb, lorb
        integer :: ioflag_, nline_
        real(kind=8) :: prec_
!f2py intent(in) a
!f2py intent(in) b
!f2py intent(in) c
!f2py intent(in) d
!f2py intent(in) e
!f2py intent(in) nline
!f2py intent(in) path
!f2py intent(in) ioflag
!f2py intent(in) prec
!f2py intent(in) mat
!f2py depend(a,b,c,d,e) mat


! initial
        ioflag_ = 1
        nline_  = 0
        prec_   = 1.0d-16
        if(present(ioflag)) ioflag_ = ioflag
        if(present(nline))  nline_  = nline
        if(present(prec))   prec_   = prec
        
! key part
        if(ioflag_ .eq. 2)then
            open(40,file = path, access="append")
        elseif(ioflag_ .eq. 1)then
            open(40,file = path, status="replace")
        else
            open(40,file = path, status="replace")
        endif
        if(nline_ .gt. 0)then
            do iorb =1, nline_
                write(40,*) ''
            enddo
        else
            continue
        endif

       do iatm= 1,a
           do iorb = 1,b
               do jorb = 1,c
                   do korb = 1,d
                       do lorb = 1,e
                           if(abs(mat(iatm,iorb,jorb,korb,lorb)) .gt. prec_)then
                               write(40,'(5i10,2(7X,f35.16))')iatm,iorb,jorb,&
                                   korb,lorb,real(mat(iatm,iorb,jorb,korb,lorb)),&
                                   aimag(mat(iatm,iorb,jorb,korb,lorb))
                           endif
                       enddo
                   enddo
               enddo
           enddo
       enddo
     close(40)


  end subroutine dump_5dc
  subroutine dump_3dc(a,b,c,mat,path,nline,ioflag,prec)
      implicit none
! the dimension of 3d matrix
     integer,    intent(in)  :: a
                 
     integer,    intent(in)  :: b
                 
     integer,    intent(in)  :: c


! the dump matrix
     complex(kind=8), intent(in) :: mat(a,b,c)

! the dump path
     character(len = *),  intent(in) :: path

! first nline lines which are introduction not data
        integer,    optional    :: nline

! io flag: 1:replace, 2:append
        integer,    optional    :: ioflag

! control prec of output
        real(kind=8),   optional    :: prec




! local variables
     integer :: iatm, ibase1, ibase2, iorb, jorb, korb
     integer :: ioflag_, nline_
        real(kind=8) :: prec_
!f2py intent(in) a
!f2py intent(in) b
!f2py intent(in) c
!f2py intent(in) nline
!f2py intent(in) path
!f2py intent(in) ioflag
!f2py intent(in) prec
!f2py intent(in) mat
!f2py depend(a,b,c) mat


! initial
        ioflag_ = 1
        nline_  = 0
        prec_   = 1.0d-16
        if(present(ioflag)) ioflag_ = ioflag
        if(present(nline))  nline_  = nline
        if(present(prec))   prec_   = prec
        
! key part
        if(ioflag_ .eq. 2)then
            open(40,file = path, access="append")
        elseif(ioflag_ .eq. 1)then
            open(40,file = path, status="replace")
        else
            open(40,file = path, status="replace")
        endif
        if(nline_ .gt. 0)then
            do iorb =1, nline_
                write(40,*) ''
            enddo
        else
            continue
        endif

       do iatm= 1,a
           do iorb = 1,b
               do jorb = 1,c
                   if(abs(mat(iatm,iorb,jorb)) .gt. prec_)then
                       write(40,'(3i10,2(7X,f35.16))')iatm,iorb,jorb,real(mat(iatm,iorb,jorb)),aimag(mat(iatm,iorb,jorb))
                   endif
               enddo
           enddo
       enddo
     close(40)
  end subroutine dump_3dc

!  
  subroutine dump_3dr(a,b,c,mat,path,nline,ioflag,prec)
      implicit none
! the dimension of 3d matrix
     integer,    intent(in)  :: a
                 
     integer,    intent(in)  :: b
                 
     integer,    intent(in)  :: c


! the dump matrix
     real(kind=8), intent(in) :: mat(a,b,c)

! the dump path
     character(len = *),  intent(in) :: path

! first nline lines which are introduction not data
        integer,    optional    :: nline

! io flag: 1:replace, 2:append
        integer,    optional    :: ioflag

! control prec of output
        real(kind=8),   optional    :: prec



! local variables
     integer :: iatm, ibase1, ibase2, iorb, jorb, korb
     integer :: ioflag_, nline_
        real(kind=8) :: prec_
!f2py intent(in) a
!f2py intent(in) b
!f2py intent(in) c
!f2py intent(in) nline
!f2py intent(in) path
!f2py intent(in) ioflag
!f2py intent(in) prec
!f2py intent(in) mat
!f2py depend(a,b,c) mat



! initial
        ioflag_ = 1
        nline_  = 0
        prec_   = 1.0d-16
        if(present(ioflag)) ioflag_ = ioflag
        if(present(nline))  nline_  = nline
        if(present(prec))   prec_   = prec
        
! key part
        if(ioflag_ .eq. 2)then
            open(40,file = path, access="append")
        elseif(ioflag_ .eq. 1)then
            open(40,file = path, status="replace")
        else
            open(40,file = path, status="replace")
        endif
        if(nline_ .gt. 0)then
            do iorb =1, nline_
                write(40,*) ''
            enddo
        else
            continue
        endif


       do iatm= 1,a
           do iorb = 1,b
               do jorb = 1,c
                   if(abs(mat(iatm,iorb,jorb)) .gt. prec_)then
                       write(40,'(3i10,5x,f35.16)')iatm,iorb,jorb,mat(iatm,iorb,jorb)
                   endif
               enddo
           enddo
       enddo
     close(40)
  end subroutine dump_3dr

  subroutine dump_2dc(a,b,mat,path,nline,ioflag,prec)
      implicit none
! the dimension of 3d matrix
     integer,    intent(in)  :: a
                 
     integer,    intent(in)  :: b
                 
! the dump matrix
     complex(kind=8), intent(in) :: mat(a,b)

! the dump path
     character(len = *),  intent(in) :: path

! first nline lines which are introduction not data
        integer,    optional    :: nline

! io flag: 1:replace, 2:append
        integer,    optional    :: ioflag

! control prec of output
        real(kind=8),   optional    :: prec

! local variables
     integer :: iatm, ibase1, ibase2, iorb, jorb, korb
     integer :: ioflag_, nline_
        real(kind=8) :: prec_
!f2py intent(in) a
!f2py intent(in) b
!f2py intent(in) nline
!f2py intent(in) path
!f2py intent(in) ioflag
!f2py intent(in) prec
!f2py intent(in) mat
!f2py depend(a,b) mat



! initial
        ioflag_ = 1
        nline_  = 0
        prec_   = 1.0d-16
        if(present(ioflag)) ioflag_ = ioflag
        if(present(nline))  nline_  = nline
        if(present(prec))   prec_   = prec
        
! key part
        if(ioflag_ .eq. 2)then
            open(40,file = path, access="append")
        elseif(ioflag_ .eq. 1)then
            open(40,file = path, status="replace")
        else
            open(40,file = path, status="replace")
        endif
        if(nline_ .gt. 0)then
            do iorb =1, nline_
                write(40,*) ''
            enddo
        else
            continue
        endif

       do iatm= 1,a
           do iorb = 1,b
               if(abs(mat(iorb,iatm)) .gt. prec_)then
                   write(40,'(2i10,2(7X,f35.16))')iorb,iatm,real(mat(iorb,iatm)),aimag(mat(iorb,iatm))
               endif
           enddo
       enddo
     close(40)
  end subroutine dump_2dc

  subroutine dump_2di(a,b,mat,path,nline,ioflag,prec)
      implicit none
! the dimension of 3d matrix
     integer,    intent(in)  :: a
                 
     integer,    intent(in)  :: b
                 
! the dump matrix
     integer, intent(in) :: mat(a,b)

! the dump path
     character(len = *),  intent(in) :: path

! first nline lines which are introduction not data
        integer,    optional    :: nline

! io flag: 1:replace, 2:append
        integer,    optional    :: ioflag

! control prec of output
        real(kind=8),   optional    :: prec


! local variables
     integer :: iatm, ibase1, ibase2, iorb, jorb, korb
     integer :: ioflag_, nline_
        real(kind=8) :: prec_
!f2py intent(in) a
!f2py intent(in) b
!f2py intent(in) nline
!f2py intent(in) path
!f2py intent(in) ioflag
!f2py intent(in) prec
!f2py intent(in) mat
!f2py depend(a,b) mat


! initial
        ioflag_ = 1
        nline_  = 0
        prec_   = 1.0d-16
        if(present(ioflag)) ioflag_ = ioflag
        if(present(nline))  nline_  = nline
        if(present(prec))   prec_   = prec
        
! key part
        if(ioflag_ .eq. 2)then
            open(40,file = path, access="append")
        elseif(ioflag_ .eq. 1)then
            open(40,file = path, status="replace")
        else
            open(40,file = path, status="replace")
        endif
        if(nline_ .gt. 0)then
            do iorb =1, nline_
                write(40,*) ''
            enddo
        else
            continue
        endif

       do iatm= 1,a
           do iorb = 1,b
               if(abs(mat(iatm,iorb)) .gt. prec_)then
                   write(40,'(2i10,7x,i10)')iatm,iorb,mat(iatm,iorb)
               endif
           enddo
       enddo
     close(40)


  end subroutine dump_2di

  subroutine dump_2dr(a,b,mat,path,nline,ioflag,prec)
      implicit none
! the dimension of 3d matrix
     integer,    intent(in)  :: a
                 
     integer,    intent(in)  :: b
                 
! the dump matrix
     real(kind=8), intent(in) :: mat(a,b)

! the dump path
     character(len = *),  intent(in) :: path

! first nline lines which are introduction not data
        integer,    optional    :: nline

! io flag: 1:replace, 2:append
        integer,    optional    :: ioflag

! control prec of output
        real(kind=8),   optional    :: prec


! local variables
     integer :: iatm, ibase1, ibase2, iorb, jorb, korb
     integer :: ioflag_, nline_
        real(kind=8) :: prec_
!f2py intent(in) a
!f2py intent(in) b
!f2py intent(in) nline
!f2py intent(in) path
!f2py intent(in) ioflag
!f2py intent(in) prec
!f2py intent(in) mat
!f2py depend(a,b) mat


! initial
        ioflag_ = 1
        nline_  = 0
        prec_   = 1.0d-16
        if(present(ioflag)) ioflag_ = ioflag
        if(present(nline))  nline_  = nline
        if(present(prec))   prec_   = prec
        
! key part
        if(ioflag_ .eq. 2)then
            open(40,file = path, access="append")
        elseif(ioflag_ .eq. 1)then
            open(40,file = path, status="replace")
        else
            open(40,file = path, status="replace")
        endif
        if(nline_ .gt. 0)then
            do iorb =1, nline_
                write(40,*) ''
            enddo
        else
            continue
        endif

       do iatm= 1,a
           do iorb = 1,b
               if(abs(mat(iatm,iorb)) .gt. prec_)then
                   write(40,'(2i10,7x,f35.16)')iatm,iorb,mat(iatm,iorb)
               endif
           enddo
       enddo
     close(40)


  end subroutine dump_2dr
  subroutine dump_1dc(a,mat,path,nline,ioflag,prec)
      implicit none
! the dimension of 3d matrix
     integer,    intent(in)  :: a
                 
! the dump matrix
     complex(kind=8), intent(in) :: mat(a)

! the dump path
     character(len = *),  intent(in) :: path

! first nline lines which are introduction not data
        integer,    optional    :: nline

! io flag: 1:replace, 2:append
        integer,    optional    :: ioflag

! control prec of output
        real(kind=8),   optional    :: prec

! local variables
     integer :: iatm, ibase1, ibase2, iorb, jorb, korb
     integer :: ioflag_, nline_
        real(kind=8) :: prec_
!f2py intent(in) a
!f2py intent(in) nline
!f2py intent(in) path
!f2py intent(in) ioflag
!f2py intent(in) prec
!f2py intent(in) mat
!f2py depend(a) mat


! initial
        ioflag_ = 1
        nline_  = 0
        prec_   = 1.0d-16
        if(present(ioflag)) ioflag_ = ioflag
        if(present(nline))  nline_  = nline
        if(present(prec))   prec_   = prec
        
! key part
        if(ioflag_ .eq. 2)then
            open(40,file = path, access="append")
        elseif(ioflag_ .eq. 1)then
            open(40,file = path, status="replace")
        else
            open(40,file = path, status="replace")
        endif
        if(nline_ .gt. 0)then
            do iorb =1, nline_
                write(40,*) ''
            enddo
        else
            continue
        endif

       do iatm= 1,a
               if(abs(mat(iatm)) .gt. prec_)then
                   write(40,'(i10,2(7x,f35.16))')iatm,real(mat(iatm)),aimag(mat(iatm))
               endif
       enddo
     close(40)
  end subroutine dump_1dc

  subroutine dump_1dr(a,mat,path,nline,ioflag,prec)
      implicit none
! the dimension of 3d matrix
     integer,    intent(in)  :: a
                 
! the dump matrix
     real(kind=8), intent(in) :: mat(a)

! the dump path
     character(len = *),  intent(in) :: path

! first nline lines which are introduction not data
        integer,    optional    :: nline

! io flag: 1:replace, 2:append
        integer,    optional    :: ioflag

! control prec of output
        real(kind=8),   optional    :: prec


! local variables
     integer :: iatm, ibase1, ibase2, iorb, jorb, korb
     integer :: ioflag_, nline_
        real(kind=8) :: prec_
!f2py intent(in) a
!f2py intent(in) nline
!f2py intent(in) path
!f2py intent(in) ioflag
!f2py intent(in) prec
!f2py intent(in) mat
!f2py depend(a) mat


! initial
        ioflag_ = 1
        nline_  = 0
        prec_   = 1.0d-16
        if(present(ioflag)) ioflag_ = ioflag
        if(present(nline))  nline_  = nline
        if(present(prec))   prec_   = prec
        
! key part
        if(ioflag_ .eq. 2)then
            open(40,file = path, access="append")
        elseif(ioflag_ .eq. 1)then
            open(40,file = path, status="replace")
        else
            open(40,file = path, status="replace")
        endif
        if(nline_ .gt. 0)then
            do iorb =1, nline_
                write(40,*) ''
            enddo
        else
            continue
        endif

       do iatm= 1,a
           if(abs(mat(iatm)) .gt. prec_)then
               write(40,'(i10,7x,f35.16)')iatm,mat(iatm)
           endif
       enddo
     close(40)


  end subroutine dump_1dr
  subroutine dump_1di(a,mat,path,nline,ioflag,prec)
      implicit none
! the dimension of 3d matrix
     integer,    intent(in)  :: a
                 
! the dump matrix
     integer,    intent(in) :: mat(a)

! the dump path
     character(len = *),  intent(in) :: path

! first nline lines which are introduction not data
        integer,    optional    :: nline

! io flag: 1:replace, 2:append
        integer,    optional    :: ioflag

! control prec of output
        real(kind=8),   optional    :: prec


! local variables
     integer :: iatm, ibase1, ibase2, iorb, jorb, korb
     integer :: ioflag_, nline_
        real(kind=8) :: prec_
!f2py intent(in) a
!f2py intent(in) nline
!f2py intent(in) path
!f2py intent(in) ioflag
!f2py intent(in) prec
!f2py intent(in) mat
!f2py depend(a) mat


! initial
        ioflag_ = 1
        nline_  = 0
        prec_   = 1.0d-16
        if(present(ioflag)) ioflag_ = ioflag
        if(present(nline))  nline_  = nline
        if(present(prec))   prec_   = prec
        
! key part
        if(ioflag_ .eq. 2)then
            open(40,file = path, access="append")
        elseif(ioflag_ .eq. 1)then
            open(40,file = path, status="replace")
        else
            open(40,file = path, status="replace")
        endif
        if(nline_ .gt. 0)then
            do iorb =1, nline_
                write(40,*) ''
            enddo
        else
            continue
        endif

       do iatm= 1,a
           if(abs(mat(iatm)) .gt. prec_)then
               write(40,'(i10,7x,I10)')iatm,mat(iatm)
           endif
       enddo
     close(40)


  end subroutine dump_1di
  subroutine dump_0dar(aa,mat,path,nline,ioflag)
      implicit none
                 
! the dump matrix
     character(len = *),  intent(in) :: aa

     real(8),             intent(in) :: mat

! the dump path
     character(len = *),  intent(in) :: path

! first nline lines which are introduction not data
        integer,    optional    :: nline

! io flag: 1:replace, 2:append
        integer,    optional    :: ioflag


! local variables
     integer :: iatm, ibase1, ibase2, iorb, jorb, korb
     integer :: ioflag_, nline_
!f2py intent(in) aa
!f2py intent(in) nline
!f2py intent(in) path
!f2py intent(in) ioflag
!f2py intent(in) mat


! initial
        ioflag_ = 1
        nline_  = 0
        if(present(ioflag)) ioflag_ = ioflag
        if(present(nline))  nline_  = nline
        
! key part
        if(ioflag_ .eq. 2)then
            open(40,file = path, access="append")
        elseif(ioflag_ .eq. 1)then
            open(40,file = path, status="replace")
        else
            open(40,file = path, status="replace")
        endif
        if(nline_ .gt. 0)then
            do iorb =1, nline_
                write(40,*) ''
            enddo
        else
            continue
        endif

        write(40,'(a,4x,f35.16)') aa, mat
        close(40)


  end subroutine dump_0dar
  subroutine dump_0da(mat,path,nline,ioflag)
      implicit none
                 
! the dump matrix
     character(len = *),  intent(in) :: mat

! the dump path
     character(len = *),  intent(in) :: path

! first nline lines which are introduction not data
        integer,    optional    :: nline

! io flag: 1:replace, 2:append
        integer,    optional    :: ioflag


! local variables
     integer :: iatm, ibase1, ibase2, iorb, jorb, korb
     integer :: ioflag_, nline_
!f2py intent(in) nline
!f2py intent(in) path
!f2py intent(in) ioflag
!f2py intent(in) mat

! initial
        ioflag_ = 1
        nline_  = 0
        if(present(ioflag)) ioflag_ = ioflag
        if(present(nline))  nline_  = nline
        
! key part
        if(ioflag_ .eq. 2)then
            open(40,file = path, access="append")
        elseif(ioflag_ .eq. 1)then
            open(40,file = path, status="replace")
        else
            open(40,file = path, status="replace")
        endif
        if(nline_ .gt. 0)then
            do iorb =1, nline_
                write(40,*) ''
            enddo
        else
            continue
        endif

        write(40,'(a)') mat
        close(40)


  end subroutine dump_0da
  subroutine dump_0di(mat,path,nline,ioflag)
      implicit none
                 
! the dump matrix
     integer, intent(in) :: mat

! the dump path
     character(len = *),  intent(in) :: path

! first nline lines which are introduction not data
        integer,    optional    :: nline

! io flag: 1:replace, 2:append
        integer,    optional    :: ioflag


! local variables
     integer :: iatm, ibase1, ibase2, iorb, jorb, korb
     integer :: ioflag_, nline_
!f2py intent(in) nline
!f2py intent(in) path
!f2py intent(in) ioflag
!f2py intent(in) mat

! initial
        ioflag_ = 1
        nline_  = 0
        if(present(ioflag)) ioflag_ = ioflag
        if(present(nline))  nline_  = nline
        
! key part
        if(ioflag_ .eq. 2)then
            open(40,file = path, access="append")
        elseif(ioflag_ .eq. 1)then
            open(40,file = path, status="replace")
        else
            open(40,file = path, status="replace")
        endif
        if(nline_ .gt. 0)then
            do iorb =1, nline_
                write(40,*) ''
            enddo
        else
            continue
        endif

        write(40,'(I8)') mat
        close(40)


  end subroutine dump_0di

  subroutine dump_0dr(mat,path,nline,ioflag)
      implicit none
                 
! the dump matrix
     real(kind=8), intent(in) :: mat

! the dump path
     character(len = *),  intent(in) :: path

! first nline lines which are introduction not data
        integer,    optional    :: nline

! io flag: 1:replace, 2:append
        integer,    optional    :: ioflag


! local variables
     integer :: iatm, ibase1, ibase2, iorb, jorb, korb
     integer :: ioflag_, nline_
!f2py intent(in) nline
!f2py intent(in) path
!f2py intent(in) ioflag
!f2py intent(in) mat

! initial
        ioflag_ = 1
        nline_  = 0
        if(present(ioflag)) ioflag_ = ioflag
        if(present(nline))  nline_  = nline
        
! key part
        if(ioflag_ .eq. 2)then
            open(40,file = path, access="append")
        elseif(ioflag_ .eq. 1)then
            open(40,file = path, status="replace")
        else
            open(40,file = path, status="replace")
        endif
        if(nline_ .gt. 0)then
            do iorb =1, nline_
                write(40,*) ''
            enddo
        else
            continue
        endif

        write(40,'(f35.16)') mat
        close(40)


  end subroutine dump_0dr
  subroutine dump_0dc(mat,path,nline,ioflag)
      implicit none
! the dump matrix
     complex(kind=8), intent(in) :: mat

! the dump path
     character(len = *),  intent(in) :: path

! first nline lines which are introduction not data
        integer,    optional    :: nline

! io flag: 1:replace, 2:append
        integer,    optional    :: ioflag

! local variables
     integer :: iatm, ibase1, ibase2, iorb, jorb, korb
     integer :: ioflag_, nline_
!f2py intent(in) nline
!f2py intent(in) path
!f2py intent(in) ioflag
!f2py intent(in) mat

! initial
        ioflag_ = 1
        nline_  = 0
        if(present(ioflag)) ioflag_ = ioflag
        if(present(nline))  nline_  = nline
        
! key part
        if(ioflag_ .eq. 2)then
            open(40,file = path, access="append")
        elseif(ioflag_ .eq. 1)then
            open(40,file = path, status="replace")
        else
            open(40,file = path, status="replace")
        endif
        if(nline_ .gt. 0)then
            do iorb =1, nline_
                write(40,*) ''
            enddo
        else
            continue
        endif

        write(40,'(f35.16)')real(mat),aimag(mat)
     close(40)
  end subroutine dump_0dc

