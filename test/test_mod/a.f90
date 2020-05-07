!program main
!    integer :: a1, a2
!    a1 = 1
!    call test1(a1,a2)
!end program
subroutine test1(a,b)
    implicit none
    integer,  intent(in)  ::  a
    integer,  intent(out) ::  b
    !
    b = a + 1
    print *, 'a=',a
    print *, 'b=',b
end subroutine
