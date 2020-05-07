!
  !
! read vpms in many body space
  subroutine rtgw_read_vpmatc(num_vpm, num_config, phimat)

      implicit none

! number of vpms
      integer,  intent(in)     ::   num_vpm

! number of orbitals
      integer,  intent(in)     ::   num_config

! fmat
      complex(kind=8), intent(out)     ::   phimat(num_vpm, num_config, num_config)

! local variables
      integer    :: iwan, jwan, ivpm, kk, ii, jj
      integer    :: iost
      real(kind=8)   :: tmp1, tmp2

! logical 
      logical    :: exists,fish

!f2py   intent(in)    num_vpm
!f2py   intent(in)    num_config
!f2py   intent(out)   phimat
!f2py   depend(num_vpm, num_config)  phimat

! key part
  phimat = dcmplx(0.d0,0.d0)
  iost = 0
  fish = .true.
  inquire(file = 'Gutz_4.vpm.nphi', exist = exists)
  if(exists)then
      open(24,file='Gutz_4.vpm.nphi',status = 'old')
      read(24,*)
      read(24,*)
      read(24,*)
      do while(fish)
          read(24,*,iostat = iost)kk, ii, jj, tmp1, tmp2
          if(iost .lt. 0)then
              fish = .false.
          elseif(iost .eq. 0)then
!              print *,kk,ii,jj,tmp1,tmp2
              phimat(kk, ii, jj) = dcmplx(tmp1, tmp2)
          else
              close(24)
              stop "read Gutz_4.vpm.nphi wrong!"
          endif
      enddo
      close(24)
  else
      stop "NO INPUT(Gutz_4.vpm.nphi)"
  endif
!      inquire(file = 'Gutz_4.vpm.nphi', exist = exists)
!      if(exists)then
!          open(24,file='Gutz_4.vpm.nphi',status='old')
!          read(24,*)
!          read(24,*)
!          read(24,*)
!          do ivpm = 1,num_vpm
!              do iwan = 1,num_config
!                  do jwan = 1,num_config
!                      read(24,*)kk, ii, jj, tmp1, tmp2
!                      if((kk .ne. ivpm) .or. (ii .ne. jwan) .or. (jj .ne. iwan))then
!                          stop "sth wrong when read Gutz_4.vpm.nphi in rtgw_HB_solver"
!                      endif
!                      phimat(ivpm, jwan, iwan) = dcmplx(tmp1,tmp2)
!                  enddo
!              enddo
!          enddo
!          close(24)
!      else
!          stop "NO INPUT(Gutz_4.vpm.nphi) FOR rtgw_HB_solver SUBROUTAIN"
!      endif

  end subroutine rtgw_read_vpmatc



! destination : read two dimension complex vector form filepath
! @sypeng
subroutine read_2dc(nline, xdm, ydm, filepath, Amat)
    implicit none

! first nline lines which are introduction not data
    integer,            intent(in)  ::   nline

! the xdim of the two-dimension vector
    integer,            intent(in)  ::   xdm

! the ydim of the two-dimension vector
    integer,            intent(in)  ::   ydm

! filepath or name of the input data
    character(len=*),   intent(in)  ::   filepath

! output : 1-D real vector
    complex(kind=8),    intent(out) ::   Amat(xdm, ydm)


! auxiliary variables
    integer       :: i, j, iost
    real(kind=8)  :: tmp1, tmp2
    logical       :: exists, fish

!f2py   intent(in)     nline
!f2py   intent(in)     xdm
!f2py   intent(in)     ydm
!f2py   intent(in)     filepath
!f2py   intent(out)    Amat
!f2py   depend(xdm,ydm)  Amat

! initial 
    exists = .false.
    fish   = .true.
    iost   = 0
    Amat   = dcmplx(0.d0, 0.d0)
    
! key part
    inquire(file = filepath, exist = exists )
    if(exists)then
        open(24, file = filepath, status = "old")
        if(nline .gt. 0)then
            do i = 1, nline
                read(24, *)
            enddo
        elseif(nline .eq. 0)then
            continue
        else
            print *,"[nline is ignored] ",nline," is minus integer, how many lines do you want to ignore in the front of ", filepath
        endif
        do while(fish)
            read(24,*,iostat = iost)i, j, tmp1, tmp2
            if(iost .lt. 0)then
                fish = .false.
            elseif(iost .eq. 0)then
                if(i .le. xdm .and. j .le. ydm)then
                    Amat(i, j) = dcmplx(tmp1, tmp2)
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

end subroutine read_2dc


! destination : read one dimension real vector form filepath
! @sypeng
subroutine read_1dr(nline, xdm, filepath, Amat)
    implicit none

! first nline lines which are introduction not data
    integer,            intent(in)  ::   nline

! the length of the one-dimension vector
    integer,            intent(in)  ::   xdm

! filepath or name of the input data
    character(len=*),   intent(in)  ::   filepath

! output : 1-D real vector
    real(kind=8),       intent(out) ::   Amat(xdm)


! auxiliary variables
    integer       :: i, iost
    real(kind=8)  :: tmp
    logical       :: exists, fish

!f2py   intent(in)     nline
!f2py   intent(in)     xdm
!f2py   intent(in)     filepath
!f2py   intent(out)    Amat
!f2py   depend(xdm)    Amat

! initial 
    exists = .false.
    fish   = .true.
    iost   = 0
    Amat   = 0.d0
    
! key part
    inquire(file = filepath, exist = exists )
    if(exists)then
        open(24, file = filepath, status = "old")
        if(nline .gt. 0)then
            do i = 1, nline
                read(24, *)
            enddo
        elseif(nline .eq. 0)then
            continue
        else
            print *,"[nline is ignored] ",nline," is minus integer, how many lines do you want to ignore in the front of ", filepath
        endif
        do while(fish)
            read(24,*,iostat = iost)i, tmp
            if(iost .lt. 0)then
                fish = .false.
            elseif(iost .eq. 0)then
                if(i .le. xdm)then
                    Amat(i) = tmp
                else
                    print *,"the index:",i,"in ",filepath," exceed the boundry of input dimension:",xdm
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
subroutine read_1dc(nline, xdm, filepath, Amat)
    implicit none

! first nline lines which are introduction not data
    integer,            intent(in)  ::   nline

! the length of the one-dimension vector
    integer,            intent(in)  ::   xdm

! filepath or name of the input data
    character(len=*),   intent(in)  ::   filepath

! output : 1-D real vector
    complex(kind=8),    intent(out) ::   Amat(xdm)


! auxiliary variables
    integer       :: i, iost
    real(kind=8)  :: tmp1, tmp2
    logical       :: exists, fish

!f2py   intent(in)     nline
!f2py   intent(in)     xdm
!f2py   intent(in)     filepath
!f2py   intent(out)    Amat
!f2py   depend(xdm)    Amat

! initial 
    exists = .false.
    fish   = .true.
    iost   = 0
    Amat   = dcmplx(0.d0, 0.d0)
    
! key part
    inquire(file = filepath, exist = exists )
    if(exists)then
        open(24, file = filepath, status = "old")
        if(nline .gt. 0)then
            do i = 1, nline
                read(24, *)
            enddo
        elseif(nline .eq. 0)then
            continue
        else
            print *,"[nline is ignored] ",nline," is minus integer, how many lines do you want to ignore in the front of ", filepath
        endif
        do while(fish)
            read(24,*,iostat = iost)i, tmp1, tmp2
            if(iost .lt. 0)then
                fish = .false.
            elseif(iost .eq. 0)then
                if(i .le. xdm)then
                    Amat(i) = dcmplx(tmp1, tmp2)
                else
                    print *,"the index:",i,"in ",filepath," exceed the boundry of input dimension:",xdm
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
subroutine read_2dr(nline, xdm, ydm, filepath, Amat)
    implicit none

! first nline lines which are introduction not data
    integer,            intent(in)  ::   nline

! the xdim of the two-dimension vector
    integer,            intent(in)  ::   xdm

! the ydim of the two-dimension vector
    integer,            intent(in)  ::   ydm

! filepath or name of the input data
    character(len=*),   intent(in)  ::   filepath

! output : 2-D real vector
    real(kind=8),       intent(out) ::   Amat(xdm, ydm)


! auxiliary variables
    integer       :: i, j, iost
    real(kind=8)  :: tmp
    logical       :: exists, fish

!f2py   intent(in)     nline
!f2py   intent(in)     xdm
!f2py   intent(in)     ydm
!f2py   intent(in)     filepath
!f2py   intent(out)    Amat
!f2py   depend(xdm,ydm)    Amat

! initial 
    exists = .false.
    fish   = .true.
    iost   = 0
    Amat   = 0.d0
    
! key part
    inquire(file = filepath, exist = exists )
    if(exists)then
        open(24, file = filepath, status = "old")
        if(nline .gt. 0)then
            do i = 1, nline
                read(24, *)
            enddo
        elseif(nline .eq. 0)then
            continue
        else
            print *,"[nline is ignored] ",nline," is minus integer, how many lines do you want to ignore in the front of ", filepath
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



! destination : read one dimension complex vector form filepath
! @sypeng
subroutine read_2dcr(nline, xdm, ydm, filepath, Amat)
    implicit none

! first nline lines which are introduction not data
    integer,            intent(in)  ::   nline

! the xdim of the two-dimension vector
    integer,            intent(in)  ::   xdm

! the ydim of the two-dimension vector
    integer,            intent(in)  ::   ydm

! filepath or name of the input data
    character(len=*),   intent(in)  ::   filepath

! output : 1-D real vector
    complex(kind=8),    intent(out) ::   Amat(xdm, ydm)


! auxiliary variables
    integer       :: i, j, iost
    real(kind=8)  :: tmp1, tmp2
    logical       :: exists, fish

!f2py   intent(in)     nline
!f2py   intent(in)     xdm
!f2py   intent(in)     ydm
!f2py   intent(in)     filepath
!f2py   intent(out)    Amat
!f2py   depend(xdm,ydm)    Amat

! initial 
    exists = .false.
    fish   = .true.
    iost   = 0
    Amat   = dcmplx(0.d0, 0.d0)
    
! key part
    inquire(file = filepath, exist = exists )
    if(exists)then
        open(24, file = filepath, status = "old")
        if(nline .gt. 0)then
            do i = 1, nline
                read(24, *)
            enddo
        elseif(nline .eq. 0)then
            continue
        else
            print *,"[nline is ignored] ",nline," is minus integer, how many lines do you want to ignore in the front of ", filepath
        endif
        do while(fish)
            read(24,*,iostat = iost)i, j, tmp1
            if(iost .lt. 0)then
                fish = .false.
            elseif(iost .eq. 0)then
                if(i .le. xdm .and. j .le. ydm)then
                    Amat(i, j) = dcmplx(tmp1, 0.d0)
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

end subroutine read_2dcr



! destination : read three dimension real vector form filepath
! @sypeng
subroutine read_3dr(nline, xdm, ydm, zdm, filepath, Amat)
    implicit none

! first nline lines which are introduction not data
    integer,            intent(in)  ::   nline

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
    integer       :: i, j, k, iost
    real(kind=8)  :: tmp
    logical       :: exists, fish

!f2py   intent(in)     nline
!f2py   intent(in)     xdm
!f2py   intent(in)     ydm
!f2py   intent(in)     zdm
!f2py   intent(in)     filepath
!f2py   intent(out)    Amat
!f2py   depend(xdm,ydm,zdm)    Amat

! initial 
    exists = .false.
    fish   = .true.
    iost   = 0
    Amat   = 0.d0
    
! key part
    inquire(file = filepath, exist = exists )
    if(exists)then
        open(24, file = filepath, status = "old")
        if(nline .gt. 0)then
            do i = 1, nline
                read(24, *)
            enddo
        elseif(nline .eq. 0)then
            continue
        else
            print *,"[nline is ignored] ",nline," is minus integer, how many lines do you want to ignore in the front of ", filepath
        endif
        do while(fish)
            read(24,*,iostat = iost)i, j, k, tmp
            if(iost .lt. 0)then
                fish = .false.
            elseif(iost .eq. 0)then
                if(i .le. xdm .and. j .le. ydm .and. k .le. zdm)then
                    Amat(i, j, k) = tmp
                else
                    print *,"the index:",i,"or",j,"or",k,"in ",filepath,&
                        " exceed the boundry of input dimension:",xdm,"or",ydm,"or",zdm
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
subroutine read_3dc(nline, xdm, ydm, zdm, filepath, Amat)
    implicit none

! first nline lines which are introduction not data
    integer,            intent(in)  ::   nline

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
    integer       :: i, j, k, iost
    real(kind=8)  :: tmp1, tmp2
    logical       :: exists, fish

!f2py   intent(in)     nline
!f2py   intent(in)     xdm
!f2py   intent(in)     ydm
!f2py   intent(in)     zdm
!f2py   intent(in)     filepath
!f2py   intent(out)    Amat
!f2py   depend(xdm,ydm,zdm)    Amat

! initial 
    exists = .false.
    fish   = .true.
    iost   = 0
    Amat   = dcmplx(0.d0, 0.d0)
    
! key part
    inquire(file = filepath, exist = exists )
    if(exists)then
        open(24, file = filepath, status = "old")
        if(nline .gt. 0)then
            do i = 1, nline
                read(24, *)
            enddo
        elseif(nline .eq. 0)then
            continue
        else
            print *,"[nline is ignored] ",nline," is minus integer, how many lines do you want to ignore in the front of ", filepath
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
                        " exceed the boundry of input dimension:",xdm,"or",ydm,"or",zdm
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

