!  FirstOrderUpwind.f90 
!
!  FUNCTIONS:
!  MAIN - Entry point of console application.
!  Author      : Xingguang Zhou
!  Organization: Xi`an Jiaotong University - NuTHeL
!  Date        : 2020/10/7
!

!****************************************************************************
!
!  PROGRAM: FirstOrderUpwind
!
!  PURPOSE:  Calculate the first order upwind scheme on -----
!            1-D steady state without internal heat source problem.
!
!****************************************************************************

    ! generate the grid, and setup some parameters.
    module SETUP
    implicit none
        private dx !The distance of each grid
        real*8 :: dx
        !define the Node structure.
        type :: Node
            real*8 :: x         !location
            real*8 :: temperature
            real*8 :: lambda    !heat conduction coefficient
            real*8 :: density
            real*8 :: velocity
        end type
        !define the Interface structure.
        type :: Surface
            real*8 :: x         !location
            real*8 :: lambda
            real*8 :: density
            real*8 :: velocity
        end type
    contains
        !begin to generate the grid with Practice B 
        subroutine Grid(NodeGroup, SurfaceGroup, N, length)
        implicit none
            type(Node), intent(inout)    :: NodeGroup(:)
            type(Surface), intent(inout) :: SurfaceGroup(:)
            real*8, intent(in)           :: length
            integer, intent(in)          :: N            !the number of the grids
            integer, parameter           :: fileid = 10
            integer                      :: i
            character(len=20)            :: filename='location.txt'
            logical                      :: alive
            dx = length / N
            !begin to generate the surface location
            SurfaceGroup(1)%x = 0
            do i=1, N
                SurfaceGroup(i+1)%x = SurfaceGroup(i)%x + dx
            end do
            !begin to generate the node location
            NodeGroup(1)%x = SurfaceGroup(1)%x
            NodeGroup(N+2)%x = SurfaceGroup(N+1)%x
            do i=1, N
                NodeGroup(i+1)%x = (SurfaceGroup(i)%x + SurfaceGroup(i+1)%x) / 2.D0
            end do
            !check the file status
            inquire(file=filename, exist=alive)
            if(alive) then
                write(*, *) "The Grid message file has generated, updating the data now..."
            else
                write(*, *) "The Grid message file has not benn generated, generating the data file now..."
            end if
            !begin to write the file
            open(unit=fileid, file=filename)
            write(fileid, *) "Node Loca:     Surface Loca:"
            do i=1, N+1
                write(fileid, "(F8.4, 8X, F8.4)") NodeGroup(i)%x, SurfaceGroup(i)%x
            end do
            write(fileid, "(F8.4, 8X, F8.4)") NodeGroup(N+2)%x
            close(fileid)
            return 
        end subroutine
        !
        !begin to generate the physical parameters to the Node and Surface.
        !
        subroutine Property(NodeGroup, SurfaceGroup)
        implicit none
            type(Node), intent(inout)    :: NodeGroup(:)
            type(Surface), intent(inout) :: SurfaceGroup(:)
            integer                      :: N               !the number of node
            integer                      :: i
            real*8                       :: const_lambda    !the constant heat conduciton parameter
            real*8                       :: const_density   !the constant density
            real*8                       :: const_velocity  !the constant velocity
            logical                      :: IsConstant
            logical                      :: alive
            integer, parameter           :: fileid = 10
            character(len=20)            :: filename='property.txt'
            N = size(NodeGroup, 1)
            write(*, *) "Constant Property?"
            read(*, *) IsConstant
            if(IsConstant == .false.) then
                !TODO: .....
                write(*, *) "Please note the property!"
                read(*, *)
            end if
            ! begin to read the constant property
            write(*, *) "Please input the heat conduction coefficient:"
            read(*, *) const_lambda
            write(*, *) "Please input the fluid velocity:"
            read(*, *) const_velocity
            write(*, *) "Please input the fluid density:"
            read(*, *) const_density
            !Only care about the parameters of the INTERNAL NODE!
            do i=2, N-1
                NodeGroup(i)%lambda = const_lambda
                NodeGroup(i)%velocity = const_velocity
                NodeGroup(i)%density = const_density
            end do
            !
            !all surface property should be given!
            !
            do i=1, N-1
                SurfaceGroup(i)%lambda = const_lambda
                SurfaceGroup(i)%velocity = const_velocity
                SurfaceGroup(i)%density = const_density
            end do
            !check the file status
            inquire(file=filename, exist=alive)
            if(alive) then
                write(*, *) "Property file is ready, updating the data now..."
            else
                write(*, *) "Property file is not ready, writing the data now..."
            end if
            !begin to write the file
            !we should pay more attention to the surface group property.
            open(unit=fileid, file=filename)
            write(fileid, *) "Node HC:       Surface HC:"
            do i=1, N-1
                write(fileid, "(F8.4, 8X, F8.4)") NodeGroup(i)%lambda, SurfaceGroup(i)%lambda
            end do
            write(fileid, *) "Node density:       Surface density:"
            do i=1, N-1
                write(fileid, "(F8.4, 8X, F8.4)") NodeGroup(i)%density, SurfaceGroup(i)%density
            end do
            write(fileid, *) "Node V:       Surface V:"
            do i=1, N-1
                write(fileid, "(F8.4, 8X, F8.4)") NodeGroup(i)%velocity, SurfaceGroup(i)%velocity
            end do
            close(fileid)
            return
        end subroutine
        ! TDMA coefficient matrix generation.
        !
        ![ap1  ae1                                   ]
        !|aw2  ap2  ae2                              |
        !|     aw3  ap3  ae3                         |
        !|          aw4  ap4  ae4                    |
        !|                ...                        |
        !|                  ...                      |
        !|                    ...                    | 
        !|                      awn-1  apn-1  aen-1  |
        ![                               awn  apn    ]
        !
        !Generate the coefficient matrix of TDMA
        !because of non-source, then there is no 'b' inside the matrix function.
        subroutine GenerateModulus(AE, AP, AW, SurfaceGroup)
        implicit none
            real*8, intent(inout)     :: AE(:)
            real*8, intent(inout)     :: AP(:)
            real*8, intent(inout)     :: AW(:)
            type(Surface), intent(in) :: SurfaceGroup(:)
            integer                   :: N
            integer                   :: i
            logical                   :: alive
            character(len=20)         :: filename='COEFF.txt'
            integer, parameter        :: fileid=10
            N = size(SurfaceGroup,1)
            !because of the fluid velocity has two direction of 1-D axis ( v>0 or v<0 ),
            !so we need to distinguish the two style upwind scheme.
            if(SurfaceGroup(2)%velocity > 0) then
                goto 1000
            else
                goto 2000
            end if
            !<1>. if v>0:
1000        AW(1) = 0.D0         !triangle-diagonal matrix condition
            do i=2, N+1
                AW(i) = SurfaceGroup(i-1)%lambda/dx + SurfaceGroup(i-1)%density*SurfaceGroup(i-1)%velocity
            end do
            AE(N+1) = 0.D0          !triangle-diagonal matrix condition
            do i=1, N
                AE(i) = SurfaceGroup(i)%lambda/dx
            end do
            goto 3000
            !<2>. if v<0
            !begin to calculate the array AE
2000        AW(1) = 0.D0
            do i=2, N+1
                AW(i) = SurfaceGroup(i-1)%lambda/dx
            end do
            !reverse the array AW
            !AW = AW(N+1:1:-1)
            !begin to calculate the array AE
            AE(N+1) = 0.D0
            do i=1, N
                AE(i) = SurfaceGroup(i)%lambda/dx + SurfaceGroup(i)%density*dabs(SurfaceGroup(i)%velocity)
            end do
            !reverse the array AE
            !AE = AE(N+1:1:-1)
            !Get array AP, AP = AW + AE, use the Fortran built-in array operation.
3000        AP = AW + AE
            !check the file status
            inquire(file=filename, exist=alive)
            if(alive) then
                write(*, *) "The COEFF file is ready, updating the data now..."
            else
                write(*, *) "The COEFF file is not ready, writing the data now..."
            end if
            !begin to write the COEFF file
            open(unit=fileid, file=filename)
            write(fileid, *) "COEFF OF TDMA:"
            write(fileid, *) "  AW               AP                AE"
            do i=1, N+1
                write(fileid, "(F8.4,8X,F8.4,8X,F8.4)") AW(i), AP(i), AE(i)
            end do
            close(fileid)
            return
        end subroutine
        !print some message of the node and interface(surface).
        subroutine PrintMessage(NodeGroup, SurfaceGroup)
        implicit none
            type(Node), intent(in)    :: NodeGroup(:)
            type(Surface), intent(in) :: SurfaceGroup(:)
            integer                   :: i, N
            N = size(SurfaceGroup, 1)
            write(*, *) "*********调用打印信息函数*********"
            write(*, *)
            write(*, *) "节点坐标：      界面坐标："
            do i=1, N
                write(*, "(F8.4,8X,F8.4)") NodeGroup(i)%x, SurfaceGroup(i)%x
            end do
            write(*, "(F8.4)") NodeGroup(N+1)%x
            write(*, *) "*********************************"
            !TODO: print more message
            return
        end subroutine
    end module
    !****************************************************************
    !TDMA solution module
    !****************************************************************
    module TDMA
    use SETUP
    implicit none
    contains
        subroutine SolvePQ(P, Q, A, B, C, t1, tm1, v_direction)
        implicit none
            real*8, intent(inout) :: P(:)
            real*8, intent(inout) :: Q(:)
            real*8, intent(in)    :: A(:)
            real*8, intent(in)    :: B(:)
            real*8, intent(in)    :: C(:)
            real*8, intent(in)    :: t1
            real*8, intent(in)    :: tm1
            logical, intent(in)   :: v_direction !record the fluid flow direction.
            real*8                :: coeff=0.D0  !record the temp parameter.
            integer               :: fileid=10
            character(len=20)     :: filename='PQ.txt'
            integer               :: i
            integer               :: N
            N = size(A, 1)                       !此时的A的元素个数与节点个数相同
            if(v_direction) then ! V > 0
                !The first kind boundary condition:
                P(1) = 0
                Q(1) = t1
            else                 ! V < 0
                P(1) = 0
                Q(1) = t1
            end if
            do i=2, N
                coeff = A(i)-C(i)*P(i-1)
                P(i) = B(i) / coeff
                Q(i) = (C(i)*Q(i-1)) / coeff
            end do
            !begin to write the file.
            open(unit=fileid, file=filename)
            write(fileid, *) "P         Q"
            do i=1, N
                write(fileid, "(F8.4, 8X, F8.4)") P(i), Q(i)
            end do
            close(fileid)
            return
        end subroutine
        !Solve the temperature of each node.
        subroutine SolveT(P, Q, NodeGroup)
        implicit none
            real*8, intent(inout)     :: P(:)              !P和Q的元素个数现在与节点数据相同，但是求解的时候用不到这么多的PQ
            real*8, intent(inout)     :: Q(:)
            type(Node), intent(inout) :: NodeGroup(:)
            integer                   :: N
            integer                   :: i
            N = size(NodeGroup, 1)
            !从倒数第二个节点开始回代求解温度
            do i=N-1, 2, -1
                NodeGroup(i)%temperature = P(i)*NodeGroup(i+1)%temperature + Q(i)
            end do
            return
        end subroutine
    end module
    !
    !The entry of the program
    !
    program main
    use SETUP
    use TDMA
    implicit none
        real*8                     :: length !the length of this physical problem.
        real*8, allocatable        :: A(:)
        real*8, allocatable        :: B(:)
        real*8, allocatable        :: C(:)
        real*8, allocatable        :: P(:)
        real*8, allocatable        :: Q(:)
        type(Node), allocatable    :: NodeGroup(:)
        type(Surface), allocatable :: SurfaceGroup(:)
        integer                    :: i
        integer                    :: N
        integer, parameter         :: fileid = 10
        logical                    :: alive
        character(len=20)          :: filename = 'temperature.csv'
        logical                    :: v_direction
        write(*, *) "Please input the length of this problem:"
        read(*, *) length
        write(*, *) "How many volume-domain you want to generate?"
        read(*, *) N
        !allocate the memory to the allocatable arrays.
        allocate(A(N+2))
        allocate(B(N+2))
        allocate(C(N+2))
        allocate(P(N+2))
        allocate(Q(N+2))
        allocate(NodeGroup(N+2))
        allocate(SurfaceGroup(N+1))
        !begin to generator the grid
        call Grid(NodeGroup, SurfaceGroup, N, length)
        !print the grid message
        call PrintMessage(NodeGroup, SurfaceGroup)
        !set the first kind boundary condition
        write(*, *) "Please input the left node temperature:"
        read(*, *) NodeGroup(1)%temperature
        write(*, *) "Please input the right node temperature:"
        read(*, *) NodeGroup(N+2)%temperature
        !begin to set the property.
        call Property(NodeGroup, SurfaceGroup)
        call GenerateModulus(C, A, B, SurfaceGroup)
        P(:) = 0.D0
        Q(:) = 0.D0
        write(*, *) "Fluid flow direction is on the x-axis positive direction?"
        read(*, *) v_direction
        !begin to solve the matrix P and Q
        call SolvePQ(P, Q, A, B, C, NodeGroup(1)%temperature, NodeGroup(N+2)%temperature, v_direction)
        !begin to solve the T
        call SolveT(P, Q, NodeGroup)
        write(*, *) "Temperature:", (NodeGroup(i)%temperature, i=1, N+2)
        !check the file status
        inquire(file=filename, exist=alive)
        if(alive) then
            write(*, *) "The temperature data file is ready, updating the data now...."
        else
            write(*, *) "The temperature data file not ready, writing the data now...."
        end if
        !begin to write the data
        open(unit=fileid, file=filename)
        do i=1, N+2
            write(fileid, *) NodeGroup(i)%temperature
        end do
        close(fileid)
        !release the memory which has benn allocated.
        deallocate(A)
        deallocate(B)
        deallocate(C)
        deallocate(P)
        deallocate(Q)
        deallocate(NodeGroup)
        deallocate(SurfaceGroup)
    end program
    