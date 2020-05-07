
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
!f2py intent(in) ioflag
!f2py intent(in) prec
!f2py intent(in) path
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

! destination : read two dimension real vector form filepath
! @sypeng
subroutine read_2dr(xdm, ydm, Amat, filepath, nline)
    implicit none

! the xdim of the two-dimension vector
    integer,            intent(in)  ::   xdm

! the ydim of the two-dimension vector
    integer,            intent(in)  ::   ydm

! output : 2-D real vector
    real(kind=8),       intent(out) ::   Amat(xdm, ydm)

! filepath or name of the input data
    character(len=*),   intent(in)  ::   filepath

! first nline lines which are introduction not data
    integer,            optional    ::   nline


! auxiliary variables
    integer       :: i, j, iost, nline_
    real(kind=8)  :: tmp
    logical       :: exists, fish

!f2py intent(in) xdm
!f2py intent(in) ydm
!f2py intent(in) nline
!f2py intent(in) filepath
!f2py intent(out) Amat
!f2py depend(xdm,ydm) Amat


! initial 
    exists = .false.
    fish   = .true.
    iost   = 0
    Amat   = 0.d0
    nline_ = 0

    if(present(nline)) nline_ = nline
    
    
! key part
    inquire(file = filepath, exist = exists )
    if(exists)then
        open(24, file = filepath, status = "old")
        if(nline_ .gt. 0)then
            do i = 1, nline_
                read(24, *)
            enddo
        elseif(nline_ .eq. 0)then
            continue
        else
           print *,"[nline is ignored] ",nline_," is minus integer, how many lines do you want to ignore in the front of ", filepath
        endif
        do while(fish)
            read(24,*,iostat = iost)i, j, tmp
            if(iost .lt. 0)then
                fish = .false.
            elseif(iost .eq. 0)then
                if(i .le. xdm .and. j .le. ydm)then
                    Amat(i, j) = tmp
                else
                    print *,"the index:",i,"or",j,"in ",filepath," exceed the boundry of input dimension:",xdm,"or",ydm
                    stop
                endif
            else
                close(24)
                print *,"something wrong when reading",filepath
                stop
            endif
        enddo
        close(24)
    else
        print *, filepath,"doesn't exist!"
        stop
    endif

end subroutine read_2dr

