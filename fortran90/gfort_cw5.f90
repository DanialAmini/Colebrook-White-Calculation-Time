

!cygwin terminal, run as admin
!cd "C:\Users\beny\Downloads\Colebrook-White-Calculation-Time-master\Colebrook-White-Calculation-Time-master"
!cd C:\gfortran_files
!gfortran gfort_cw4.f90 
!--> this wont' work, you have to compile with openmp because of 
!openmp identifiers in the program
!gfortran -O3 -mavx -fopenmp gfort_cw4.F90
!./a.exe


! the best one by a long margin is manual index multithreading
!but maybe not worth the trouble?
!manual hand-tuned code for multithreading


!second best is simple loops
!the easiest to code too


!2d matrix: each cols = problem size/cores, number of rows=threads
!vectorization as fast as loops if manual multithreading is done


!Dlog is log for real(8), real(8) has 8 bytes = double precision (float64)
!dlog and log have the same speed


!it doesn't even ork I think
!openMP fortran gide: https://curc.readthedocs.io/en/latest/programming/OpenMP-Fortran.html


!http://fortranwiki.org/fortran/show/Real+precision
!integer, parameter :: dp=kind(0.d0) !I think this gets the double size?
!integer, parameter :: dp = selected_real_kind(15, 307)    !fortran 90 also
!use, intrinsic :: iso_fortran_env !2008 notation
!integer, parameter :: dp = REAL64 !2008 notation
!integer, parameter :: dp = kind(1.d0)     !fortran 90
!integer, parameter :: dp = kind(0.d0)     !fortran 90

!for some reason "allocate" should be done 
!after all of the variables are defined

!for some reason it's using single precision like NOOOO

program cw_func

    use omp_lib
    

    

    implicit none
    integer :: thread_id

integer, parameter :: dp = kind(0.d0)     !fortran 90

    INTEGER :: cr, cm,c1,c2,n, i,t ,tnr,j, nn_1,n8
    REAL(kind=dp) :: rate,r

    real(kind=dp), allocatable :: Re(:),epsD(:)
    real(kind=dp), allocatable :: z_(:),z_hat(:) , y1(:), y2(:), A(:)
    real(kind=dp), allocatable :: Re_(:,:),epsD_(:,:) , A_(:,:)
    
    real(kind=dp) :: k1,k2,k3,k4,k5,k6,k7


    n=10
    n=10000000
    
    allocate(A(n))
    allocate(Re(n))
    allocate(epsD(n))
    allocate(z_(n))
    allocate(z_hat(n))
    allocate(y1(n))
    allocate(y2(n))

    allocate(Re_(n8,nn_1))
    allocate(epsD_(n8,nn_1))
    allocate(A_(n8,nn_1))


    CALL system_clock(count_rate=cr)
    CALL system_clock(count_max=cm)  
    rate = REAL(cr)

    call omp_set_num_threads(8)
    t = omp_get_num_threads()


    k1=4000
    k2=1e8
    k3=1e-6
    k4=0.05
    call srand(86456)
    do i = 1, n
        call random_number(r)
        Re(i) = k1*dexp(r*dlog(k2/k1))
        call random_number(r)
        epsD(i) = k3*(k4/k3)**r
    end do

    Re(1)=156423.1564812345561
    epsD(1)=0.0011564564538731584


    CALL SYSTEM_CLOCK(c1)
    !$omp parallel do
    do i = 1, n
        z_(i)=func_sj(Re(i),epsD(i))
    end do
    !$omp end parallel do
    CALL SYSTEM_CLOCK(c2)
    print "('#2.00000 system_clock:', f15.10,'s')", (c2-c1)/rate


    k1=3.7
    k2=5.74
    k3=0.9
    k4=10
    k5=0.25
    CALL SYSTEM_CLOCK(c1)
    !$omp parallel do
    do i = 1, n
        A(i)=LOG(epsD(i)/k1+k2/Re(i)**k3)/log(k4)
        A(i)=k5/(A(i)*A(i))
    end do
    !$omp end parallel do
    CALL SYSTEM_CLOCK(c2)
    print "('#2 system_clock:', f15.10,'s')", (c2-c1)/rate
    !second best


    write(*,*) 'Re,epsD,f'
    do i=1,10
          write(*,*) Re(i),',',epsD(i),',',A(i)
    end do



    contains

        function func_sj(Re__,epsD__) result(z__)
            !use iso_fortran_env, only: int64, dp !2008

            real(kind=dp) :: Re__, epsD__, z__,k1,k2,k3,k4,k5
            integer :: i
            
            k1=10
            k2=3.7
            k3=5.74
            k4=0.9
            k5=0.25
            
            z__=LOG(epsD__/k2+k3/Re__**k4)/log(k1)
            z__=k5/(z__*z__)

        end function func_sj
        





    end program cw_func