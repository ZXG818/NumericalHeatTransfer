    !作者：周星光
    !组织：Xi`an Jiaotong University - NuTHeL
    !日期：2020/9/30
    !程序：数值传热学中一维离散方程的三对角矩阵求解算法-TDMA算法实现。
    !注意：依据我的理解，TDMA算法由于不同节点相关系数不知道，只知道边界条件与相关的递推关系式，
    !      所以我们即使写出来三对角矩阵，则非零元素也是未知的，故不能直接使用矩阵的相关行变换来
    !      来进行矩阵消元，依据陶文铨教授上课讲的递推关系式，我们可以利用初始条件来逐个逐个的节点
    !      来进行递推，最后求得所有未知项的系数，此时可以进行高斯消元或者进一步利用递推关系式来求得各个节点的温度T。
    !修改：2020年10月5日，对一维稳态导热问题的网格的划分进行了更正，
    !      使其满足区域离散方法B的网格要求。
    !
    !****************************************************************
    !一维稳态导热问题初始条件设置模块（区域离散方法B，均分网格）
    !常物性方程
    !****************************************************************
    !
    module SETUP
    implicit none
        !网格间距dx，指界面与界面的距离
        private delta_x !设置为私有变量
        real*8 :: delta_x
        !尝试增加一维稳态网格划分，使用笛卡尔坐标系，网格为均匀网格，方法为有限体积方法-FVM(Finite Volume Method)
        !定义一个节点结构体
        type :: Node
            real*8 :: x             !节点坐标
            real*8 :: temperature   !节点温度
            real*8 :: lambda        !节点处的导热系数
        end type
        !定义一个界面结构体
        type :: Surface             !由于interface是Fortran中的关键字，所以采用surface来表示界面
            real*8 :: x             !界面坐标
            real*8 :: lambda        !界面处的当量导热系数
        end type
    contains
        !开始对一维网格进行划分
        !2020/10/5，对一维网格划分进行修正，使其符合区域离散方法B的网格要求。
        subroutine Grid(NodeGroup, SurfaceGroup, N, length)
        implicit none
            type(Node), intent(inout)    :: NodeGroup(:)              !用于存储节点的数组
            type(Surface), intent(inout) :: SurfaceGroup(:)           !用于存储界面的数组
            integer, intent(in)          :: N                         !要进行划分的控制体积的数目
            real*8, intent(in)           :: length                    !问题边界的整体长度
            integer                      :: i                         
            integer, parameter           :: fileid = 10               !文件ID
            character(len=20)            :: filename='location.txt'   !文件名
            logical                      :: alive
            delta_x = length / N
            !网格区域数目，界面数目，节点数目之间的关系是依次递增一
            !接下来对网格进行划分，记录节点坐标与界面坐标
            !先对界面进行划分
            SurfaceGroup(1)%x = 0
            do i=1, N
                SurfaceGroup(i+1)%x = SurfaceGroup(i)%x + delta_x
            end do
            !界面坐标划分完毕后，再对节点坐标进行划分，
            !在区域离散方法B中，节点除了问题边界处在和界面坐标相同外，其余均在两个界面的中心位置
            NodeGroup(1)%x = SurfaceGroup(1)%x
            NodeGroup(N+2)%x = SurfaceGroup(N+1)%x
            do i=1, N
                NodeGroup(i+1)%x = (SurfaceGroup(i)%x+SurfaceGroup(i+1)%x) / 2.D0
            end do
            !检查文件状态
            inquire(file=filename, exist=alive)
            if(alive) then
                write(*, *) "物性状态文件已经存在，现重新更新文件内容...."
            else
                write(*, *) "物性状态文件不存在，现开始创建文件...."
            end if
            !开始写入文件
            open(unit=fileid, file=filename)
            write(fileid, *) "节点坐标：    界面坐标："
            do i=1, N+1
                write(fileid, "(F8.4, 8X, F8.4)") NodeGroup(i)%x, SurfaceGroup(i)%x
            end do
            write(fileid, "(F8.4, 8X, F8.4)") NodeGroup(N+2)%x
            close(fileid)
            return
        end subroutine
        !开始对网格和节点处的物性进行赋值与求解
        subroutine Property(NodeGroup, SurfaceGroup)
        implicit none
            type(Node), intent(inout)    :: NodeGroup(:)   
            type(Surface), intent(inout) :: SurfaceGroup(:)
            integer                      :: N          !用于记录节点的个数
            integer                      :: i
            real*8                       :: const_lambda              !常物性的导热系数
            logical                      :: IsConstant !用于判断用户的导热系数是否为常数
            logical                      :: alive
            integer, parameter           :: fileid=10  !文件ID
            character(len=20)            :: filename='property.txt'
            N = size(NodeGroup, 1)       !节点个数
            write(*, *)"该问题中，导热系数是否为常数？"
            read(*, *) IsConstant
            if(IsConstant == .true.) then
                goto 1001
            end if
            !
            !由于TDMA在第一类边界条件下是从第二个节点到倒数第二个节点之间进行运算的，
            !所以只需要对NodeGroup的第二个节点与倒数第二个节点之间进行赋物性值即可
            !同样，处于问题边界处的界面也不用考虑赋物性值。
            !因为只有这样，才可以使用形成三对角矩阵来使用TDMA算法进行求解。
            !
            !让用户开始输入节点处的导热系数
            write(*, *)"请输入节点处的导热系数:"
            do i=1, N
                write(*, *) NodeGroup(i)%lambda
            end do
            !开始进行求解界面处的当量导热系数，
            !使用调和平均法来求解界面处的导热系数
            do i=1, N-1
                SurfaceGroup(i)%lambda = 2*NodeGroup(i)%lambda*NodeGroup(i+1)%lambda/ &
                                         (NodeGroup(i)%lambda+NodeGroup(i+1)%lambda)
            end do
            goto 1002
1001        write(*, *)"请输入常物性的导热系数："
            read(*, *) const_lambda
            do i=1, N-1
                NodeGroup(i)%lambda = const_lambda
                SurfaceGroup(i)%lambda = const_lambda
            end do
            NodeGroup(N)%lambda = const_lambda
            !给倒数第二个节点赋导热系数
            NodeGroup(N-1)%lambda = const_lambda
            !检查文件状态
1002        inquire(file=filename, exist=alive)
            if(alive) then
                write(*, *) "物性状态文件已经存在，现重新更新文件内容...."
            else
                write(*, *) "物性状态文件不存在，现开始创建文件...."
            end if
            !开始写入文件
            open(unit=fileid, file=filename)
            write(fileid, *) "节点导热率：  界面导热率："
            do i=1, N-1
                write(fileid, "(F8.4,8X,F8.4)") NodeGroup(i)%lambda, SurfaceGroup(i)%lambda
            end do
            write(fileid, "(F8.4)") NodeGroup(N)%lambda
            close(fileid)
            return 
        end subroutine
        !
        !TODO:还应当查找导热系数与温度的关系式，然后通过迭代来确定lambda（瞬态问题下，较复杂，需要迭代）
        !
        !仍然需要注意的是：TDMA的消元与迭代过程依旧是在第2个节点到倒数第2个节点之间完成的，
        !因为只有这样，一维稳态方程的矩阵才可以是三对角矩阵
        !同时，第1个节点与最后1个节点的作用是提供第一类边界条件，来生成第二个节点到倒数第二个节点的系数,
        !由此来组成系数矩阵。
        !
        !=============================================================================
        !此时，界面和节点的坐标以及其物性参数均已设置完毕,                              
        !接下来需要依据离散完成的公式来对各个节点上的温度进行求解。                     
        !=============================================================================
        !通过界面信息生成统一形式中aP*TP = aE*TE + aW*TW + b 的aP、aE、aW以及b
        !其中A包含了AE、AW，aP用AP表示
        !生成的三对角矩阵示意图：
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
        !注意：aw1=0，aen=0，为边界条件，故在三对角矩阵中没有写入，
        !同时这样的三对角矩阵也符合也符合计算方法上的三对角矩阵形式。
        !注意此时AE,AP,AW,b的元素个数应当与节点的元素个数相同，但是求解的时候则不会用到这么多的元素。
        !
        subroutine GenerateModulus(AE, AP, AW, b, SurfaceGroup, Source)
        implicit none
            real*8, intent(inout)       :: AE(:)
            real*8, intent(inout)       :: AW(:)
            real*8, intent(inout)       :: AP(:)
            real*8, intent(inout)       :: b(:)
            real*8, intent(in)          :: Source     !内热源项，这里为了简化问题，将内热源设定为常数
            type(Surface), intent(in)   :: SurfaceGroup(:)
            integer                     :: N          !记录数组个数
            integer                     :: i
            logical                     :: alive
            character(len=20)           :: filename='COEFF.txt'
            integer, parameter          :: fileid = 10
            !我们不需要和边界重合的界面，所以我们将N设置为size(SurfaceGroup, 1) - 1
            N = size(SurfaceGroup, 1)
            !生成AW信息
            AW(1) = 0                                 !边界条件
            do i=2, N+1
                AW(i) = SurfaceGroup(i-1)%lambda / delta_x
            end do
            !生成AE信息
            AE(N+1) = 0                                 !边界条件
            do i=1, N
                AE(i) = SurfaceGroup(i)%lambda / delta_x
            end do
            !
            !TODO:AP与源强有关，需要进一步讨论，这里为了简化问题，将内热源设定为常数
            !
            !Fortran内置的数组运算
            AP = AE + AW          !本来还要在减去SpApdx，但是默认内热源为恒定源强，所以Sp=0
            b = Source*delta_x
            !全部数据生成完毕，下面开始写入文件中
            !检查文件状态
            inquire(file=filename, exist=alive)
            if(alive) then
                write(*, *)"统一方程各项系数文件已经存在，正在重新写入数据..."
            else
                write(*, *)"统一方程各项数据文件不存在，正在创建..."
            end if
            open(unit=fileid, file=filename)
            write(fileid, *) "lambda:     AW           AP               AE                 b"
            do i=1, N+1
                write(fileid, "(12X,F8.4,8X,F8.4,8X,F8.4,8X,F8.4)") AW(i), AP(i), AE(i), b(i)
            end do
            close(fileid)
            return
        end subroutine
        !打印节点和界面的一些信息
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
            !TODO:打印更多的信息
            return
        end subroutine
    end module
    !****************************************************************
    !TDMA求解模块
    !****************************************************************
    module TDMA
    use SETUP
    implicit none
    contains
        !先求解P和Q，求出P和Q以后，我们就可以利用高斯消元解矩阵或者利用迭代来求解T。
        !给TDMA方法添加第一类边界条件，即给出第一个节点的温度t1和最后一个节点的温度tm1
        subroutine SolvePQ(P, Q, A, B, C, D, t1, tm1)
        implicit none
            real*8, intent(inout) :: P(:)
            real*8, intent(inout) :: Q(:)
            real*8, intent(in)    :: A(:)           !Ti前系数
            real*8, intent(in)    :: B(:)           !Ti+1前系数
            real*8, intent(in)    :: C(:)           !Ti-1前系数
            real*8, intent(in)    :: D(:)
            real*8, intent(in)    :: t1
            real*8, intent(in)    :: tm1
            real*8                :: coeff = 0.D0   !记录中间系数
            integer               :: i
            integer               :: fileid=10
            character(len=20)     :: filename='PQ.txt'
            integer               :: length         !记录A,B,C,D,P,Q数组的长度，也就是需要求解的内部节点的个数
            length = size(A, 1)
            !首先依据第一类边界条件，来求得P1与Q1
            P(1) = 0
            Q(1) = t1
            !下面开始进行迭代来进行求解P与Q
            !使用了第一类边界条件，所以消元过程从第二个节点开始
            do i=2, length
                coeff = A(i)-C(i)*P(i-1)
                P(i) = B(i) / coeff
                Q(i) = (D(i)+C(i)*Q(i-1)) / coeff
            end do
            !开始写文件
            open(unit=fileid, file=filename)
            write(fileid, *) "P         Q"
            do i=1, length
                write(fileid, "(F8.4, 8X, F8.4)") P(i), Q(i)
            end do
            close(fileid)
            return
        end subroutine
        !求出P与Q之后，接下来开始进行回代来求解各个节点的温度
        subroutine SolveT(P, Q, NodeGroup)
        implicit none
            real*8, intent(inout)        :: P(:)
            real*8, intent(inout)        :: Q(:)
            type(Node), intent(inout)    :: NodeGroup(:)
            integer                      :: length
            integer                      :: i
            length = size(NodeGroup, 1)
            !开始从后向前进行回代求解温度
            !回代过程是从倒数第二个节点到第二个节点
            !先求倒数第二个节点的温度
            NodeGroup(length-1)%temperature = Q(length-1)
            do i=length-2, 2, -1
                NodeGroup(i)%temperature = P(i)*NodeGroup(i+1)%temperature + Q(i)
            end do
            return
        end subroutine
    end module
    !开始主函数
    program main
    use SETUP
    use TDMA
    implicit none
        real*8                    :: Source !源项
        real*8                    :: Length !问题的长度
        real*8, allocatable       :: A(:)
        real*8, allocatable       :: B(:)
        real*8, allocatable       :: C(:)
        real*8, allocatable       :: D(:)
        real*8, allocatable       :: P(:)
        real*8, allocatable       :: Q(:)
        type(Node), allocatable   :: NodeGroup(:)
        type(Surface),allocatable :: SurfaceGroup(:)
        integer                   :: i
        integer                   :: N    !记录申请空间的大小
        integer, parameter        :: fileid = 10
        logical                   :: alive
        character(len=20)         :: filename = 'temperature.csv'
        write(*, *) "请输入所描述问题的长度:"
        read(*, *) Length
        write(*, *) "请输入划分的控制体积的个数："
        read(*, *) N
        !开辟内存空间
        allocate(A(N+2))
        allocate(B(N+2))
        allocate(C(N+2))
        allocate(D(N+2))
        allocate(P(N+2))
        allocate(Q(N+2))
        allocate(NodeGroup(N+2))
        allocate(SurfaceGroup(N+1))
        !开始划分网格
        call Grid(NodeGroup, SurfaceGroup, N, Length)
        !打印网格信息
        call PrintMessage(NodeGroup, SurfaceGroup)
!DEBUG下的条件编译
!!DEC$ DEFINE __DEBUG__
!DEC$ IF DEFINED(__DEBUG__)
        !下面给数组赋初值：
        !边界条件：C(1)=0，B(10)=0
        write(*, *) "请输入Ti 项的系数："
        read(*, *) A(:)
        write(*, *) "请输入Ti+1项的系数:"
        read(*, *) B(:)
        write(*, *) "请输入Ti-1项的系数:"
        read(*, *) C(:)
        !假设该一维导热问题中没有内热源，则数组D为0
        !D = 0.D0
        write(*, *) "请输入常数项D："
        read(*, *) D(:)
!DEC$ ENDIF
        !设置第一类边界条件
        write(*, *)"请输入第一个节点的温度:"
        read(*, *) NodeGroup(1)%temperature
        write(*, *)"请输入最后一个节点的温度"
        read(*, *) NodeGroup(N+2)%temperature
        !开始设置物性
        call Property(NodeGroup, SurfaceGroup)
        write(*, *)"请输入源项："
        read(*, *) Source
        call GenerateModulus(C, A, B, D, SurfaceGroup, Source)
        !给P,Q数组初始值全部赋值为0
        P(:) = 0.D0
        Q(:) = 0.D0
        !开始进行系数P和Q的求解
        call SolvePQ(P, Q, A, B, C, D, NodeGroup(1)%temperature, NodeGroup(N+2)%temperature)
        !开始进行T的求解
        call SolveT(P, Q, NodeGroup)
        write(*, *)"温度的结果为:", (NodeGroup(i)%temperature, i=1, N+2)
        !检查文件状态
        inquire(file=filename, exist=alive)
        if(alive) then
            write(*, *) "温度数据文件已经存在，现重新写入数据进行更新...."
        else
            write(*, *) "温度数据文件不存在，现开始写入数据....."
        end if
        !开始写入数据
        open(unit=fileid, file=filename)
        do i=1, N+2
            write(fileid, *) NodeGroup(i)%temperature
        end do
        !关闭文件
        close(fileid)
        !释放内存
        deallocate(A)
        deallocate(B)
        deallocate(C)
        deallocate(D)
        deallocate(P)
        deallocate(Q)
        deallocate(NodeGroup)
        deallocate(SurfaceGroup)
    end program
    