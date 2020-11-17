!cd "C:\Users\beny\Downloads\Colebrook-White-Calculation-Time-master\Colebrook-White-Calculation-Time-master"
!gfortran -O3 -mavx -fopenmp gfort_trigonometric.F90
!./a.exe


program test_evaluate_functions

    use omp_lib
    implicit none

    !use timer_module
    integer, parameter :: dp=kind(0.d0)


    integer, parameter    :: num_iterations = 10000
    integer               :: n, i,t ,tnr,j, nn_1
    real(dp)              :: r, start, finish, tarray(2)
    real(dp)              :: a_min = -1500.00
    real(dp)              :: a_max =  1500.00
    real(dp)              :: h
    
    real(dp), allocatable :: vect1(:)
    real(dp), allocatable :: vect2(:,:)
    real(dp), allocatable :: A(:)
    character(len=30)     :: arg(1)
    !type(Timer)          :: elpTime

    ! Get the dimension from the command line
!    call getarg(1, arg(1))
!    read(arg(1), *) n

    INTEGER :: cr, cm
    REAL :: rate
    INTEGER :: c1,c2

  CALL system_clock(count_rate=cr)
  CALL system_clock(count_max=cm)
  
  
  rate = REAL(cr)


    call omp_set_num_threads(8)
    t = omp_get_num_threads()

    n=10000000

    call srand(86456)
    !call elpTime%startTimer()

    
    nn_1=16

    allocate(vect1(n))
    allocate(vect2(n,nn_1))

    allocate(A(n))


    call srand(86456)
    do i = 1, n
        call random_number(r)
        vect1(i) = r
        do j=1,nn_1
        call random_number(r)
            vect2(i,j) = r        
        end do
        call random_number(r)
        A(i)=r
    end do

    do j=1,10
        print "( 'hi' , f15.10 , 'hi' )", vect1(j)
    end do


    CALL SYSTEM_CLOCK(c1)
    !$omp parallel
        A(:) = DLOG(vect1(:))       
    !$omp end parallel
    CALL SYSTEM_CLOCK(c2)
    print "('#1 system_clock:', f15.10,'s')", (c2-c1)/rate


    CALL SYSTEM_CLOCK(c1)
    !$omp parallel do
    do i = 1, 8
        A((i-1)*125000:i*125000-1)=DLOG(vect1((i-1)*125000:i*125000-1))
    end do
    !$omp end parallel do
    CALL SYSTEM_CLOCK(c2)
    print "('#1.5 system_clock:', f15.10,'s')", (c2-c1)/rate
    ! the best one by a long margin


    CALL SYSTEM_CLOCK(c1)
    !$omp parallel do
    do i = 1, n
        A(i)=DLOG(vect1(i))
    end do
    !$omp end parallel do
    CALL SYSTEM_CLOCK(c2)
    print "('#2 system_clock:', f15.10,'s')", (c2-c1)/rate



    CALL SYSTEM_CLOCK(count=c1)
    !$omp parallel
    do j=1,nn_1
        vect2(:,j)=LOG(vect2(:,j))
    end do
    !$omp end parallel
    CALL SYSTEM_CLOCK(count=c2)
    print "('#3 system_clock:', f15.10,'s')", REAL(c2-c1)/rate
    


   
    CALL SYSTEM_CLOCK(c1)
    !!$omp parallel
        A(:) = LOG(vect1(:))    
    !!$omp end parallel
    CALL SYSTEM_CLOCK(c2)
    print "('#4 system_clock:', f15.10,'s')", (c2-c1)/rate


    print "( 'err=', f15.10 , '-' , f15.10 , 'ok' )", vect1(1), vect2(1,1)

    end program test_evaluate_functions