subroutine read_0dr(Amat,filepath, nline)
    implicit none

! first nline lines which are introduction not data
    integer,            optional    ::   nline

! filepath or name of the input data
    character(len=*),   intent(in)  ::   filepath

! output : 1-D real vector
    real(kind=8),       intent(out) ::   Amat


! auxiliary variables
    integer       :: i, iost, nline_
    real(kind=8)  :: tmp
    logical       :: exists, fish

!f2py intent(in) nline
!f2py intent(in) filepath
!f2py intent(out) Amat


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
            print *,"[nline is ignored] ",nline_," is minus integer,&
                how many lines do you want to ignore in the front of ",&
                filepath
        endif
        read(24,*,iostat = iost) tmp
        if(iost .lt. 0)then
            fish = .false.
        elseif(iost .eq. 0)then
                Amat = tmp
        else
            close(24)
            print *,"something wrong when reading",filepath
            stop
        endif
        close(24)
    else
        print *, filepath,"doesn't exist!"
        stop
    endif

end subroutine read_0dr



! destination : read one dimension complex vector form filepath
! @sypeng
subroutine read_0dc(Amat, filepath, nline)
    implicit none

! first nline lines which are introduction not data
    integer,            optional    ::   nline

! filepath or name of the input data
    character(len=*),   intent(in)  ::   filepath

! output : 1-D real vector
    complex(kind=8),    intent(out) ::   Amat


! auxiliary variables
    integer       :: i, iost, nline_
    real(kind=8)  :: tmp1, tmp2
    logical       :: exists, fish

!f2py intent(in) nline
!f2py intent(in) filepath
!f2py intent(out) Amat


! initial 
    exists = .false.
    fish   = .true.
    iost   = 0
    Amat   = dcmplx(0.d0, 0.d0)
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
            print *,"[nline is ignored] ",nline_," is minus integer,&
                how many lines do you want to ignore in the front of ",&
                filepath
        endif
        read(24,*,iostat = iost)tmp1, tmp2
        if(iost .lt. 0)then
            fish = .false.
        elseif(iost .eq. 0)then
            Amat = dcmplx(tmp1, tmp2)
        else
            close(24)
            print *,"something wrong when reading",filepath
            stop
        endif
        close(24)
    else
        print *, filepath,"doesn't exist!"
        stop
    endif

end subroutine read_0dc

! destination : read one dimension real vector form filepath
! @sypeng
subroutine read_1dr(xdm, Amat, filepath, nline)
    implicit none

! first nline lines which are introduction not data
    integer,            optional    ::   nline

! the length of the one-dimension vector
    integer,            intent(in)  ::   xdm

! filepath or name of the input data
    character(len=*),   intent(in)  ::   filepath

! output : 1-D real vector
    real(kind=8),       intent(out) ::   Amat(xdm)


! auxiliary variables
    integer       :: i, iost, nline_
    real(kind=8)  :: tmp
    logical       :: exists, fish


!f2py intent(in) xdm
!f2py intent(in) nline
!f2py intent(in) filepath
!f2py intent(out) Amat
!f2py depend(xdm) Amat


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
            print *,"[nline is ignored] ",nline_," is minus integer,&
                how many lines do you want to ignore in the front of ",&
                filepath
        endif
        do while(fish)
            read(24,*,iostat = iost)i, tmp
            if(iost .lt. 0)then
                fish = .false.
            elseif(iost .eq. 0)then
                if(i .le. xdm)then
                    Amat(i) = tmp
                else
                    print *,"the index:",i,"in ",filepath,&
                        " exceed the boundry of input dimension:",xdm
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

end subroutine read_1dr



! destination : read one dimension complex vector form filepath
! @sypeng
subroutine read_1dc(xdm, Amat, filepath, nline)
    implicit none

! first nline lines which are introduction not data
    integer,            optional    ::   nline

! the length of the one-dimension vector
    integer,            intent(in)  ::   xdm

! filepath or name of the input data
    character(len=*),   intent(in)  ::   filepath

! output : 1-D real vector
    complex(kind=8),    intent(out) ::   Amat(xdm)


! auxiliary variables
    integer       :: i, iost, nline_
    real(kind=8)  :: tmp1, tmp2
    logical       :: exists, fish


!f2py intent(in) xdm
!f2py intent(in) nline
!f2py intent(in) filepath
!f2py intent(out) Amat
!f2py depend(xdm) Amat

! initial 
    exists = .false.
    fish   = .true.
    iost   = 0
    Amat   = dcmplx(0.d0, 0.d0)
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
            print *,"[nline is ignored] ",nline_," is minus integer,&
                how many lines do you want to ignore in the front of ",&
                filepath
        endif
        do while(fish)
            read(24,*,iostat = iost)i, tmp1, tmp2
            if(iost .lt. 0)then
                fish = .false.
            elseif(iost .eq. 0)then
                if(i .le. xdm)then
                    Amat(i) = dcmplx(tmp1, tmp2)
                else
                    print *,"the index:",i,"in ",filepath,&
                        " exceed the boundry of input dimension:",xdm
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

end subroutine read_1dc


! destination : read two dimension real vector form filepath
! @sypeng
subroutine read_2dr(xdm, ydm, Amat, filepath, nline)
    implicit none

! first nline lines which are introduction not data
    integer,            optional    ::   nline

! the xdim of the two-dimension vector
    integer,            intent(in)  ::   xdm

! the ydim of the two-dimension vector
    integer,            intent(in)  ::   ydm

! filepath or name of the input data
    character(len=*),   intent(in)  ::   filepath

! output : 2-D real vector
    real(kind=8),       intent(out) ::   Amat(xdm, ydm)


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
            print *,"[nline is ignored] ",nline_," is minus integer,&
                how many lines do you want to ignore in the front of ",&
                filepath
        endif
        do while(fish)
            read(24,*,iostat = iost)i, j, tmp
            if(iost .lt. 0)then
                fish = .false.
            elseif(iost .eq. 0)then
                if(i .le. xdm .and. j .le. ydm)then
                    Amat(i, j) = tmp
                else
                    print *,"the index:",i,"or",j,"in ",filepath,&
                        " exceed the boundry of input dimension:",&
                        xdm,"or",ydm
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

! destination : read two dimension complex vector form filepath
! @sypeng
subroutine read_2dc(xdm, ydm, Amat, filepath, nline)
    implicit none

! first nline lines which are introduction not data
    integer,            optional    ::   nline

! the xdim of the two-dimension vector
    integer,            intent(in)  ::   xdm

! the ydim of the two-dimension vector
    integer,            intent(in)  ::   ydm

! filepath or name of the input data
    character(len=*),   intent(in)  ::   filepath

! output : 1-D real vector
    complex(kind=8),    intent(out) ::   Amat(xdm, ydm)


! auxiliary variables
    integer       :: i, j, iost, nline_
    real(kind=8)  :: tmp1, tmp2
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
    Amat   = dcmplx(0.d0, 0.d0)
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
            print *,"[nline is ignored] ",nline_," is minus integer,&
                how many lines do you want to ignore in the front of ",&
                filepath
        endif
        do while(fish)
            read(24,*,iostat = iost)i, j, tmp1, tmp2
            if(iost .lt. 0)then
                fish = .false.
            elseif(iost .eq. 0)then
                if(i .le. xdm .and. j .le. ydm)then
                    Amat(i, j) = dcmplx(tmp1, tmp2)
                else
                    print *,"the index:",i,"or",j,"in ",filepath,&
                        " exceed the boundry of input dimension:",&
                        xdm,"or",ydm
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

end subroutine read_2dc

! destination : read three dimension real vector form filepath
! @sypeng
subroutine read_3dr(xdm, ydm, zdm, Amat, filepath, nline)
    implicit none

! first nline lines which are introduction not data
    integer,            optional    ::   nline

! the xdim of the three-dimension vector
    integer,            intent(in)  ::   xdm

! the ydim of the three-dimension vector
    integer,            intent(in)  ::   ydm

! the zdim of the three-dimension vector
    integer,            intent(in)  ::   zdm

! filepath or name of the input data
    character(len=*),   intent(in)  ::   filepath

! output : 2-D real vector
    real(kind=8),       intent(out) ::   Amat(xdm, ydm, zdm)


! auxiliary variables
    integer       :: i, j, k, iost, nline_
    real(kind=8)  :: tmp
    logical       :: exists, fish

!f2py intent(in) xdm
!f2py intent(in) ydm
!f2py intent(in) zdm
!f2py intent(in) nline
!f2py intent(in) filepath
!f2py intent(out) Amat
!f2py depend(xdm,ydm,zdm) Amat

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
            print *,"[nline is ignored] ",nline_," is minus integer,&
                how many lines do you want to ignore in the front of ",&
                filepath
        endif
        do while(fish)
            read(24,*,iostat = iost)i, j, k, tmp
            if(iost .lt. 0)then
                fish = .false.
            elseif(iost .eq. 0)then
                if(i .le. xdm .and. j .le. ydm .and. k .le. zdm)then
                    Amat(i, j, k) = tmp
                else
                    print *,"the index:",i,"or",j,"or",k,"in ",&
                        filepath,&
                        " exceed the boundry of input dimension:",&
                        xdm,"or",ydm,"or",zdm
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

end subroutine read_3dr



! destination : read one dimension complex vector form filepath
! @sypeng
subroutine read_3dc(xdm, ydm, zdm, Amat, filepath, nline)
    implicit none

! first nline lines which are introduction not data
    integer,            optional    ::   nline

! the xdim of the three-dimension vector
    integer,            intent(in)  ::   xdm

! the ydim of the three-dimension vector
    integer,            intent(in)  ::   ydm

! the zdim of the three-dimension vector
    integer,            intent(in)  ::   zdm

! filepath or name of the input data
    character(len=*),   intent(in)  ::   filepath

! output : 3-D real vector
    complex(kind=8),    intent(out) ::   Amat(xdm, ydm, zdm)


! auxiliary variables
    integer       :: i, j, k, iost, nline_
    real(kind=8)  :: tmp1, tmp2
    logical       :: exists, fish


!f2py intent(in) xdm
!f2py intent(in) ydm
!f2py intent(in) zdm
!f2py intent(in) nline
!f2py intent(in) filepath
!f2py intent(out) Amat
!f2py depend(xdm,ydm,zdm) Amat

! initial 
    exists = .false.
    fish   = .true.
    iost   = 0
    Amat   = dcmplx(0.d0, 0.d0)
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
            print *,"[nline is ignored] ",nline_," is minus integer, &
                how many lines do you want to ignore in the front of ",&
                filepath
        endif
        do while(fish)
            read(24,*,iostat = iost)i, j, k, tmp1, tmp2
            if(iost .lt. 0)then
                fish = .false.
            elseif(iost .eq. 0)then
                if(i .le. xdm .and. j .le. ydm .and. k .le. zdm)then
                    Amat(int(i), int(j), int(k)) = dcmplx(tmp1, tmp2)
                else
                    print *,"the index:",i,"or",j,"or",k,"in ",filepath,&
                        " exceed the boundry of input dimension:",xdm,&
                        "or",ydm,"or",zdm
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

end subroutine read_3dc
! destination : read four dimension complex vector form filepath
! @sypeng
subroutine read_4dc(xdm, ydm, zdm, kdm, Amat, filepath, nline)
    implicit none

! first nline lines which are introduction not data
    integer,            optional    ::   nline

! the xdim of the four-dimension vector
    integer,            intent(in)  ::   xdm

! the ydim of the four-dimension vector
    integer,            intent(in)  ::   ydm

! the zdim of the four-dimension vector
    integer,            intent(in)  ::   zdm

! the fourth of the four-dimension vector
    integer,            intent(in)  ::   kdm

! filepath or name of the input data
    character(len=*),   intent(in)  ::   filepath

! output : 4-D real vector
    complex(kind=8),    intent(out) ::   Amat(xdm, ydm, zdm, kdm)


! auxiliary variables
    integer       :: i, j, k, m, iost, nline_
    real(kind=8)  :: tmp1, tmp2
    logical       :: exists, fish

!f2py intent(in) xdm
!f2py intent(in) ydm
!f2py intent(in) zdm
!f2py intent(in) kdm
!f2py intent(in) nline
!f2py intent(in) filepath
!f2py intent(out) Amat
!f2py depend(xdm,ydm,zdm,kdm) Amat

! initial 
    exists = .false.
    fish   = .true.
    iost   = 0
    Amat   = dcmplx(0.d0, 0.d0)
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
            print *,"[nline is ignored] ",nline_," is minus integer,&
                how many lines do you want to ignore in the front of ",&
                filepath
        endif
        do while(fish)
            read(24,*,iostat = iost)i, j, k, m, tmp1, tmp2
            if(iost .lt. 0)then
                fish = .false.
            elseif(iost .eq. 0)then
                if(i .le. xdm .and. j .le. ydm .and. k .le. zdm .and. m .le. kdm)then
                    Amat(int(i), int(j), int(k), int(m)) = dcmplx(tmp1, tmp2)
                else
                    print *,"the index:",i,"or",j,"or",k,"or",m,"in ",&
                        filepath," exceed the boundry of inputdimension:",&
                        xdm,"or",ydm,"or",zdm,"or",kdm
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

end subroutine read_4dc

! destination : read Five dimension complex vector form filepath
! @sypeng
subroutine read_5dc(xdm, ydm, zdm, kdm, ldm, Amat, filepath, nline)
    implicit none

! first nline lines which are introduction not data
    integer,            optional    ::   nline

! the xdim of the four-dimension vector
    integer,            intent(in)  ::   xdm

! the ydim of the four-dimension vector
    integer,            intent(in)  ::   ydm

! the zdim of the four-dimension vector
    integer,            intent(in)  ::   zdm

! the fourth of the four-dimension vector
    integer,            intent(in)  ::   kdm

! the fiveth of the four-dimension vector
    integer,            intent(in)  ::   ldm

! filepath or name of the input data
    character(len=*),   intent(in)  ::   filepath

! output : 4-D real vector
    complex(kind=8),    intent(out) ::   Amat(xdm, ydm, zdm, kdm, ldm)


! auxiliary variables
    integer       :: i, j, k, m, l, iost, nline_
    real(kind=8)  :: tmp1, tmp2
    logical       :: exists, fish


!f2py intent(in) xdm
!f2py intent(in) ydm
!f2py intent(in) zdm
!f2py intent(in) kdm
!f2py intent(in) ldm
!f2py intent(in) nline
!f2py intent(in) filepath
!f2py intent(out) Amat
!f2py depend(xdm,ydm,zdm,kdm,lmat) Amat

! initial 
    exists = .false.
    fish   = .true.
    iost   = 0
    Amat   = dcmplx(0.d0, 0.d0)
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
            print *,"[nline is ignored] ",nline_," is minus integer,&
                how many lines do you want to ignore in the front of ",&
                filepath
        endif
        do while(fish)
            read(24,*,iostat = iost)i, j, k, m, l, tmp1, tmp2
            if(iost .lt. 0)then
                fish = .false.
            elseif(iost .eq. 0)then
                if(i .le. xdm .and. j .le. ydm .and. k .le. &
                    zdm .and. m .le. kdm .and. l .le. ldm)then
                    Amat(int(i), int(j), int(k), int(m), int(l))&
                        = dcmplx(tmp1, tmp2)
                else
                    print *,"the index:",i,"or",j,"or",k,"or",&
                        m,"or",l,"in ",&
                        filepath," exceed the boundry of inputdimension:",&
                        xdm,"or",ydm,"or",zdm,"or",kdm,"or",ldm
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

end subroutine read_5dc
! destination : 
! @sypeng
subroutine read_2d1dc(xdm, ydm, Amat, filepath, nline, ioflag)
    implicit none

! first nline lines which are introduction not data
    integer,            optional    ::   nline

! the xdim of the three-dimension vector
    integer,            intent(in)  ::   xdm

! the ydim of the three-dimension vector
    integer,            intent(in)  ::   ydm

! filepath or name of the input data
    character(len=*),   intent(in)  ::   filepath

! output : 3-D real vector
    complex(kind=8),    intent(out) ::   Amat(xdm*ydm)

! ioflag : 1: old
    integer,            optional    ::   ioflag


! auxiliary variables
    integer       :: i, j, k, iost,ibase1,ibase2
    real(kind=8)  :: tmp1, tmp2
    logical       :: exists, fish
    integer       :: ioflag_, nline_

!f2py intent(in) xdm
!f2py intent(in) ydm
!f2py intent(in) nline
!f2py intent(in) ioflag
!f2py intent(in) filepath
!f2py intent(out) Amat
!f2py depend(xdm,ydm) Amat

! initial 
    exists  = .false.
    fish    = .true.
    iost    = 0
    Amat    = dcmplx(0.d0, 0.d0)
    ioflag_ = 1
    nline_  = 0

    if(present(ioflag)) ioflag_ = ioflag
    if(present(nline))  nline_  = nline
! key part
    inquire(file = filepath, exist = exists )
    if(exists)then
        if(ioflag_ == 1)then
            open(24, file = filepath, status = "old")
        else
            print *,"strange ioflag=",ioflag_,"in reading",filepath
            stop
        endif
        if(nline_ .gt. 0)then
            do i = 1, nline_
                read(24, *)
            enddo
        elseif(nline_ .eq. 0)then
            continue
        else
            print *,"[nline is ignored] ",nline_," is minus integer,&
                how many lines do you want to ignore in the front of ",&
                filepath
        endif
        do while(fish)
            read(24,*,iostat = iost)i, j, tmp1, tmp2
            if(iost .lt. 0)then
                fish = .false.
            elseif(iost .eq. 0)then
                if(i .le. xdm .and. j .le. ydm)then
                    ibase1 = (int(i) -1) * ydm
                    Amat(int(j)+ibase1) = dcmplx(tmp1, tmp2)
                else
                    print *,"the index:",i,"or",j,"or",filepath,&
                        " exceed the boundry of input dimension:",&
                        xdm,"or",ydm
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

end subroutine read_2d1dc
! destination : read 3d real matrix but output a 2d real matrix
! @sypeng
subroutine read_2d1dr(xdm, ydm, Amat, filepath, nline, ioflag)
    implicit none

! first nline lines which are introduction not data
    integer,            optional    ::   nline

! the xdim of the three-dimension vector
    integer,            intent(in)  ::   xdm

! the ydim of the three-dimension vector
    integer,            intent(in)  ::   ydm

! filepath or name of the input data
    character(len=*),   intent(in)  ::   filepath

! output : 2-D real vector
    real(kind=8),       intent(out) ::   Amat(xdm*ydm)

! ioflag : 1: old
    integer,            optional    ::   ioflag



! auxiliary variables
    integer       :: i, j, k, iost,ibase1,ibase2
    real(kind=8)  :: tmp
    logical       :: exists, fish
    integer       :: nline_, ioflag_

!f2py intent(in) xdm
!f2py intent(in) ydm
!f2py intent(in) nline
!f2py intent(in) ioflag
!f2py intent(in) filepath
!f2py intent(out) Amat
!f2py depend(xdm,ydm) Amat

! initial 
    exists  = .false.
    fish    = .true.
    iost    = 0
    Amat    = 0.d0
    ioflag_ = 1
    nline_  = 0
    
    if(present(ioflag)) ioflag_ = ioflag
    if(present(nline))  nline_  = nline
! key part
    inquire(file = filepath, exist = exists )
    if(exists)then
        if(ioflag_ == 1)then
            open(24, file = filepath, status = "old")
        else
            print *,"strange ioflag=",ioflag_,&
                "in reading",filepath
            stop
        endif
        if(nline_ .gt. 0)then
            do i = 1, nline_
                read(24, *)
            enddo
        elseif(nline_ .eq. 0)then
            continue
        else
            print *,"[nline is ignored] ",nline_," is minus integer,&
                how many lines do you want to ignore in the front of ",&
                filepath
        endif
        do while(fish)
            read(24,*,iostat = iost)i, j, tmp
            if(iost .lt. 0)then
                fish = .false.
            elseif(iost .eq. 0)then
                if(i .le. xdm .and. j .le. ydm )then
                    ibase1 = (int(i) -1) * ydm
                    Amat(int(j)+ibase1) = tmp
                else
                    print *,"the index:",i,"or",j,"or",filepath,&
                        " exceed the boundry of input dimension:",&
                        xdm,"or",ydm
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

end subroutine read_2d1dr
! destination : 
! @sypeng
subroutine read_3d2dc(xdm, ydm, zdm, Amat, filepath, nline, ioflag)
    implicit none

! first nline lines which are introduction not data
    integer,            optional    ::   nline

! the xdim of the three-dimension vector
    integer,            intent(in)  ::   xdm

! the ydim of the three-dimension vector
    integer,            intent(in)  ::   ydm

! the zdim of the three-dimension vector
    integer,            intent(in)  ::   zdm

! filepath or name of the input data
    character(len=*),   intent(in)  ::   filepath

! output : 3-D real vector
    complex(kind=8),    intent(out) ::   Amat(xdm*ydm, xdm*zdm)

! ioflag : 1: old
    integer,            optional    ::   ioflag


! auxiliary variables
    integer       :: i, j, k, iost,ibase1,ibase2
    real(kind=8)  :: tmp1, tmp2
    logical       :: exists, fish
    integer       :: ioflag_, nline_
!f2py intent(in) xdm
!f2py intent(in) ydm
!f2py intent(in) zdm
!f2py intent(in) nline
!f2py intent(in) ioflag
!f2py intent(in) filepath
!f2py intent(out) Amat
!f2py depend(xdm,ydm,zdm) Amat

! initial 
    exists  = .false.
    fish    = .true.
    iost    = 0
    Amat    = dcmplx(0.d0, 0.d0)
    ioflag_ = 1
    nline_  = 0

    if(present(ioflag)) ioflag_ = ioflag
    if(present(nline))  nline_  = nline
! key part
    inquire(file = filepath, exist = exists )
    if(exists)then
        if(ioflag_ == 1)then
            open(24, file = filepath, status = "old")
        else
            print *,"strange ioflag=",ioflag_,"in reading",filepath
            stop
        endif
        if(nline_ .gt. 0)then
            do i = 1, nline_
                read(24, *)
            enddo
        elseif(nline_ .eq. 0)then
            continue
        else
            print *,"[nline is ignored] ",nline_," is minus integer,&
                how many lines do you want to ignore in the front of ",&
                filepath
        endif
        do while(fish)
            read(24,*,iostat = iost)i, j, k, tmp1, tmp2
            if(iost .lt. 0)then
                fish = .false.
            elseif(iost .eq. 0)then
                if(i .le. xdm .and. j .le. ydm .and. k .le. zdm)then
                    ibase1 = (int(i) -1) * ydm
                    ibase2 = (int(i) -1) * zdm
                    Amat(int(j)+ibase1, int(k)+ibase2) = dcmplx(tmp1, tmp2)
                else
                    print *,"the index:",i,"or",j,"or",k,"in ",filepath,&
                        " exceed the boundry of input dimension:",&
                        xdm,"or",ydm,"or",zdm
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

end subroutine read_3d2dc
! destination : read 3d real matrix but output a 2d real matrix
! @sypeng
subroutine read_3d2dr(xdm, ydm, zdm, Amat, filepath, nline, ioflag)
    implicit none

! first nline lines which are introduction not data
    integer,            optional    ::   nline

! the xdim of the three-dimension vector
    integer,            intent(in)  ::   xdm

! the ydim of the three-dimension vector
    integer,            intent(in)  ::   ydm

! the zdim of the three-dimension vector
    integer,            intent(in)  ::   zdm

! filepath or name of the input data
    character(len=*),   intent(in)  ::   filepath

! output : 2-D real vector
    real(kind=8),       intent(out) ::   Amat(xdm*ydm, xdm*zdm)

! ioflag : 1: old
    integer,            optional    ::   ioflag



! auxiliary variables
    integer       :: i, j, k, iost,ibase1,ibase2
    real(kind=8)  :: tmp
    logical       :: exists, fish
    integer       :: nline_, ioflag_
!f2py intent(in) xdm
!f2py intent(in) ydm
!f2py intent(in) zdm
!f2py intent(in) nline
!f2py intent(in) ioflag
!f2py intent(in) filepath
!f2py intent(out) Amat
!f2py depend(xdm,ydm,zdm) Amat


! initial 
    exists  = .false.
    fish    = .true.
    iost    = 0
    Amat    = 0.d0
    ioflag_ = 1
    nline_  = 0
    
    if(present(ioflag)) ioflag_ = ioflag
    if(present(nline))  nline_  = nline
! key part
    inquire(file = filepath, exist = exists )
    if(exists)then
        if(ioflag_ == 1)then
            open(24, file = filepath, status = "old")
        else
            print *,"strange ioflag=",ioflag_,&
                "in reading",filepath
            stop
        endif
        if(nline_ .gt. 0)then
            do i = 1, nline_
                read(24, *)
            enddo
        elseif(nline_ .eq. 0)then
            continue
        else
            print *,"[nline is ignored] ",nline_," is minus integer,&
                how many lines do you want to ignore in the front of ",&
                filepath
        endif
        do while(fish)
            read(24,*,iostat = iost)i, j, k, tmp
            if(iost .lt. 0)then
                fish = .false.
            elseif(iost .eq. 0)then
                if(i .le. xdm .and. j .le. ydm .and. k .le. zdm)then
                    ibase1 = (int(i) -1) * ydm
                    ibase2 = (int(i) -1) * zdm
                    Amat(int(j)+ibase1, int(k)+ibase2) = tmp
                else
                    print *,"the index:",i,"or",j,"or",k,"in ",filepath,&
                        " exceed the boundry of input dimension:",&
                        xdm,"or",ydm,"or",zdm
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

end subroutine read_3d2dr
