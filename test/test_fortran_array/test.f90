subroutine test_array(n,a,b)
    implicit none
    integer, intent(in)  :: n 
    integer, intent(in)  :: a(0:n)
    integer, intent(out) :: b(0:n)

    integer  :: i

!f2py intent(in)  n
!f2py intent(in)  a
!f2py intent(out) b
!f2py depend(n)   a
!f2py depend(n)   b
    print * , 'a=',a
    print * , 'a(0)=',a(0)
    print * , 'a(1)=',a(1)
    print * , 'a(2)=',a(2)
    print * , 'a(3)=',a(3)
    print * , 'a(4)=',a(4)
    do i = 0, n
        b(i) = a(i)
    enddo
    print * , 'b=',b
    print * , 'b(0)=',b(0)
    print * , 'b(1)=',b(1)
    print * , 'b(2)=',b(2)
    print * , 'b(3)=',b(3)
    print * , 'b(4)=',b(4)
end subroutine test_array
