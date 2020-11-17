

!cygwin terminal, run as admin
!cd "C:\Users\beny\Downloads\Colebrook-White-Calculation-Time-master\Colebrook-White-Calculation-Time-master"
!cd C:\gfortran_files
!gfortran -fopenmp -O3 -mavx  gfort_cw6.f90

!--> this wont' work, you have to compile with openmp because of 
!openmp identifiers in the program

!./a.exe


! the best one by a long margin is manual index multithreading
!but maybe not worth the trouble?
!manual hand-tuned code for multithreading


!second best is simple loops
!the easiest to code too


!2d matrix: each cols = problem size/cores, number of rows=threads
!vectorization as fast as loops if manual multithreading is done


!log is log for real(dp)(8), real(dp)(8) has 8 bytes = real(dp) (float64)
!log and log have the same speed


!it doesn't even ork I think
!openMP fortran gide: https://curc.readthedocs.io/en/latest/programming/OpenMP-Fortran.html


!http://fortranwiki.org/fortran/show/real(dp)+precision
!integer, parameter :: dp=kind(0.d0) !I think this gets the double size?
!integer, parameter :: dp = selected_real(dp)_kind(15, 307)    !fortran 90 also
!use, intrinsic :: iso_fortran_env !2008 notation
!integer, parameter :: dp = real(dp)64 !2008 notation
!integer, parameter :: dp = kind(1.d0)     !fortran 90
!integer, parameter :: dp = kind(0.d0)     !fortran 90
!real(kind=real64)
!for some reason "allocate" should be done 
!after all of the variables are defined

!for some reason it's using single precision like NOOOO

program cw_func

    use omp_lib
    implicit none

    
    integer :: thread_id

    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)  !for some reason only this works!

    INTEGER :: cr, cm,c1,c2,n, i,t ,tnr,j, nn_1,n8
    real(dp) :: rate,r

    real(dp), allocatable :: Re(:),epsD(:)
    real(dp), allocatable :: z_(:),z_hat(:) , y1(:), y2(:), A(:)
    real(dp), allocatable :: Re_(:,:),epsD_(:,:) , A_(:,:)
    
    real(dp) :: k1,k2,k3,k4,k5,k6,k7


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
    rate = real(cr)

    call omp_set_num_threads(8)
    t = omp_get_num_threads()


    k1=4000_dp
    k2=1e8_dp
    k3=1e-6_dp
    k4=0.05_dp
    call srand(86456)
    do i = 1, n
        call random_number(r)
        Re(i) = k1*exp(r*log(k2/k1))
        call random_number(r)
        epsD(i) = k3*(k4/k3)**r
    end do

    Re(1)=156423.1564812345561_dp
    epsD(1)=0.0011564564538731584_dp


    CALL SYSTEM_CLOCK(c1)
    !$omp parallel do
    do i = 1, n
        z_hat(i)=func_sj(Re(i),epsD(i))
    end do
    !$omp end parallel do
    CALL SYSTEM_CLOCK(c2)
    print "('#sj func system_clock:', f15.10,'s')", (c2-c1)/rate


    CALL SYSTEM_CLOCK(c1)
    !$omp parallel do
    do i = 1, n
        z_(i)=func_cw30(Re(i),epsD(i))
    end do
    !$omp end parallel do
    CALL SYSTEM_CLOCK(c2)
    print "('#cw30 system_clock:', f15.10,'s')", (c2-c1)/rate


    CALL SYSTEM_CLOCK(c1)
    !$omp parallel do
    do i = 1, n
        z_(i)=func_cw5(Re(i),epsD(i))
    end do
    !$omp end parallel do
    CALL SYSTEM_CLOCK(c2)
    print "('#cw5 system_clock:', f15.10,'s')", (c2-c1)/rate




    write(*,*) 'Re,epsD,f'
    do i=1,10
          write(*,*) Re(i),',',epsD(i),',',z_(i),',',z_hat(i)
    end do



    contains

        function func_sj(Re,epsD) result(z)
            !use iso_fortran_env, only: int64, dp !2008

            real(dp) :: Re, epsD, z , k1,k2,k3
            
            k1=3.7_dp
            k2=5.74_dp
            k3=0.9_dp

            k4=1.3254745276195995026404165971485_dp   
            !k6=0.25*log(10.0_dp)**2
            
            z=LOG(epsD/k1+k2/Re**k3)
            z=k4/(z*z)

        end function func_sj
        

        function func_cw30(Re,epsD) result(z)

            real(dp) :: Re, epsD, z, k1,k2,k3,k4
            integer :: i
            
            k1=3.7_dp
            !k2=log(10.0_dp)/(-2.0_dp)*(2.51_dp)
            k2=-2.8897442917075273334425792756289
            !k3=(-2/log(10.0_dp))**2 
            k3=0.75444678804645571687984331986024
            
            k1=k1*epsD
            k2=k2/Re

            z=-8.1_dp

            do i=1,30
                !z=log(k1+k2*log(k1+k2*log(k1+k2*log(k1+k2*log(k1+k2*z)))))
                z=log(k1+k2*z)
            end do

            z=k3/(z*z)

        end function func_cw30

        function func_cw5(Re,epsD) result(z)

            real(dp) :: Re, epsD, z, k1,k2,k3,k4
            
            k1=3.7_dp
            !k2=log(10.0_dp)/(-2.0_dp)*(2.51_dp)
            k2=-2.8897442917075273334425792756289
            !k3=(-2/log(10.0_dp))**2 
            k3=0.75444678804645571687984331986024
            
            k1=k1*epsD
            k2=k2/Re

            z=-8.1_dp

            do i=1,5
                z=log(k1+k2*log(k1+k2*log(k1+k2*log(k1+k2*log(k1+k2*z)))))
            end do

            z=k3/(z*z)

        end function func_cw5




    end program cw_func