    program main
    implicit none
        integer :: a(4) = (/ 1, 2, 3, 4 /)
        integer :: b(4)
        b = a(size(a,1):1:-1)
        write(*, *) b
        write(*, *) a
        a = a(size(a,1):1:-1)
        write(*, *) a
    end program
    