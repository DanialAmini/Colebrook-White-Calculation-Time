//open cygwin terminal and change the directory
//cd "C:\Users\beny\Downloads\Colebrook-White-Calculation-Time-master\Colebrook-White-Calculation-Time-master"
// g++ -Ofast -mavx -fopenmp gcc_cw6.cpp


// g++ -Ofast -mavx -fopenmp cw_gcc_9sep2020.cpp
//  export OMP_NUM_THREADS=8

// ./a.exe
// cygwin
// (OMP_NUM_THREADS not set) GNU (gcc, g++, gfortran)	-fopenmp
// -Ofast for optimization
// mavx for avx instructions
// https://www.math.ucla.edu/~wotaoyin/windows_coding_cygwin.html
// install the 64bit one from here: https://cygwin.com/install.html
// select these packages
//under Devel: 
//autoconf, autoconf2.1, autogen, automake, all automake1.*, bison, flex, gcc-core/fortran/g, gdb, git, gperf, help2man, libtool, make, patch, patchutils, pkg-config, swig, texinfo, texinfo-tex, mingw64-x86_64-gcc-core/g++
//under Libs: 
//libgcc1, lapack, liblapack-devel/doc/0, libopenmpi/-devel/12/cxx1
//under Perl: perl
//under Utils: diffutils, dos2unix
//under Web: wget
// select all packages for g++, gcc, x86_64, don't select 32 ones, don't select 861 ones
// select all openMP, MPI, blas & lapack ones
// regular timer does't work for multithreading with clock 
// loop unrolling makes performance worse (contrary to intel c++)


// https://xyce.sandia.gov/documentation/BuildingGuide.html#ubuntuPreReq

// ubuntu is debian.
// Debian Linux and Variants
// All required packages are available in these systems' default repositories.

// gcc (see the Trilinos note)
// g++
// gfortran
// make
// cmake
// bison
// flex
// libfl-dev
// libfftw3-dev
// libsuitesparse-dev
// libblas-dev
// liblapack-dev
// libtool
// If you are building the parallel version of Xyce you also need:

// libopenmpi-dev
// openmpi-bin
// Install these packages using sudo apt-get install, or a graphical package manager. The "-dev" packages will also pull additional packages in as dependencies.




// for windows: 
// Cygwin
// gcc-core
// gcc-g++
// gcc-fortran
// make
// cmake (needed for configuring Trilinos)
// lapack, liblapack0, liblapack-devel
// libsuitesparseconfig-devel
// libamd-devel
// bison
// flex
// fftw and libfftw3-devel
// libtool
// m4
// Many of these packages will automatically install additional packages that they require.

// The Xyce team has never attempted to build the parallel version of Xyce on Windows. While Cygwin does provide a version of openmpi, it is our understanding that ParMETIS does not currently build easily on Cygwin


#include <iostream>
#include <chrono>
#include <sys/types.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
//#include <math.h>

void ReEpsInitialize(double* Re__, double* epsD__, int N, int type_);
void funct_cw_iter(double* Re__, double* epsD__, double* z__, int NN, int N_iter, int N_repeat);
double MaxAbsRelErr(double* z_, double* z_hat, int NN);
double err_subs_CW(double* Re__, double* epsD__, double* z_, int NN);
void funct_approx(double* Re__, double* epsD__, double* z, int NN, int N_repeat, int method_);
void funct_cw_newton_iter(double* Re__, double* epsD__, double* z, int NN, int N_iter, int N_repeat);
void funct_zp_cw(double* Re, double* z, double* zp, int NN, int N_iter, int N_repeat);

int main()
{
    int N;
    N = 3164 - 1;
    N = 3164 - 1;
    //N = 10;
    N = 3164 - 1; //10 mill
    N = 400 - 1; //160000
    N = 200 - 1; //40000
    N = 78 - 1; //6000
    N = 32 - 1; //1000
    N = 10 - 1; //100
	N=3164-1;
//N=1000;

    struct timespec start1, finish1;
    struct timespec start00, finish00;
double elapsed;
clock_gettime(CLOCK_MONOTONIC, &start00);






    const int NN = (N + 1) * (N + 1);
    int i, j;
    int N_iter, N_repeat;
    double err_;
    int geometric_series;  
    geometric_series = 0;//arithmetic series
    geometric_series = 1;//geometric series

    double* Re__ = new double[NN];
    double* epsD__ = new double[NN];
    double* z__ = new double[NN];
    double* z4__ = new double[NN];
    double* zp = new double[NN];

    double* z_hat = new double[NN];

    N_repeat = 1;
 

    std::cout.precision(17);

    clock_t start, finish;
    clock_t start0, finish0;
    srand(time(NULL));

    printf("matrix size=%d*%d or %d elements\n\n", N + 1, N + 1, NN);

    //initialize Re & epsD
start0=clock();
    clock_gettime(CLOCK_MONOTONIC, &start1);
    ReEpsInitialize(Re__, epsD__, N, geometric_series); 
    clock_gettime(CLOCK_MONOTONIC, &finish1);
    printf("initiliaze Re & epsD duration= %lf s\n\n",  (double)((finish1.tv_sec - start1.tv_sec)+(finish1.tv_nsec - start1.tv_nsec) / 1000000000.0)  );


    //cw 50 iter
    clock_gettime(CLOCK_MONOTONIC, &start1);       N_iter = 50;       funct_cw_iter(Re__, epsD__, z__, NN, N_iter,N_repeat);        clock_gettime(CLOCK_MONOTONIC, &finish1);
    printf("Colebrook-White, %d iterations, time= %lf s, rep=%d\n", N_iter, (double)((finish1.tv_sec - start1.tv_sec)+(finish1.tv_nsec - start1.tv_nsec) / 1000000000.0) ,N_repeat);
    printf("err substituting in CW equation = %g\n\n", err_subs_CW(Re__, epsD__, z__, NN));




    //cw 3,4,5,30 iters
    clock_gettime(CLOCK_MONOTONIC, &start1);       N_iter = 3;      funct_cw_iter(Re__, epsD__, z4__, NN, N_iter, N_repeat);       clock_gettime(CLOCK_MONOTONIC, &finish1);
    printf("Colebrook-White, %d iterations, time= %lf s, rep=%d\n", N_iter, (double)((finish1.tv_sec - start1.tv_sec)+(finish1.tv_nsec - start1.tv_nsec) / 1000000000.0)  , N_repeat);
    printf("Max Rel Abs Error=%g%%\n\n", MaxAbsRelErr(z__, z4__, NN));

    clock_gettime(CLOCK_MONOTONIC, &start1);       N_iter = 4;        funct_cw_iter(Re__, epsD__, z4__, NN, N_iter, N_repeat);        clock_gettime(CLOCK_MONOTONIC, &finish1);
    printf("Colebrook-White, %d iterations, time= %lf s, rep=%d\n", N_iter, (double)((finish1.tv_sec - start1.tv_sec)+(finish1.tv_nsec - start1.tv_nsec) / 1000000000.0) , N_repeat);
    printf("Max Rel Abs Error=%g%%\n\n", MaxAbsRelErr(z__, z4__, NN));

    clock_gettime(CLOCK_MONOTONIC, &start1);       N_iter = 5;        funct_cw_iter(Re__, epsD__, z4__, NN, N_iter, N_repeat);       clock_gettime(CLOCK_MONOTONIC, &finish1);
    printf("Colebrook-White, %d iterations, time= %lf s, rep=%d\n", N_iter, (double)((finish1.tv_sec - start1.tv_sec)+(finish1.tv_nsec - start1.tv_nsec) / 1000000000.0) , N_repeat);
    printf("Max Rel Abs Error=%g%%\n\n", MaxAbsRelErr(z__, z4__, NN));

    clock_gettime(CLOCK_MONOTONIC, &start1);      N_iter = 30;        funct_cw_iter(Re__, epsD__, z4__, NN, N_iter, N_repeat);       clock_gettime(CLOCK_MONOTONIC, &finish1);
    printf("Colebrook-White, %d iterations, time= %lf s, rep=%d\n", N_iter, (double)((finish1.tv_sec - start1.tv_sec)+(finish1.tv_nsec - start1.tv_nsec) / 1000000000.0) , N_repeat);
    printf("Max Rel Abs Error=%g%%\n\n", MaxAbsRelErr(z__, z4__, NN));



    //newton & halley

    clock_gettime(CLOCK_MONOTONIC, &start1);        N_iter = 2;      funct_cw_newton_iter(Re__, epsD__, z4__, NN, 1, N_repeat);        clock_gettime(CLOCK_MONOTONIC, &finish1);
    printf("Colebrook-White Newton, %d iterations, time= %lf s, rep=%d\n", N_iter, (double)((finish1.tv_sec - start1.tv_sec)+(finish1.tv_nsec - start1.tv_nsec) / 1000000000.0)  , N_repeat);
    printf("Max Rel Abs Error=%g%%\n\n", MaxAbsRelErr(z__, z4__, NN));

    clock_gettime(CLOCK_MONOTONIC, &start1);      N_iter = 5;       funct_cw_newton_iter(Re__, epsD__, z4__, NN, 2, N_repeat);       clock_gettime(CLOCK_MONOTONIC, &finish1);
    printf("Colebrook-White Newton, %d iterations, time= %lf s, rep=%d\n", N_iter, (double)((finish1.tv_sec - start1.tv_sec)+(finish1.tv_nsec - start1.tv_nsec) / 1000000000.0) , N_repeat);
    printf("Max Rel Abs Error=%g%%\n\n", MaxAbsRelErr(z__, z4__, NN));

    clock_gettime(CLOCK_MONOTONIC, &start1);    N_iter = 2;      funct_cw_newton_iter(Re__, epsD__, z4__, NN, 3, N_repeat);        clock_gettime(CLOCK_MONOTONIC, &finish1);
    printf("Colebrook-White Hallay, %d iterations, time= %lf s, rep=%d\n", N_iter, (double)((finish1.tv_sec - start1.tv_sec)+(finish1.tv_nsec - start1.tv_nsec) / 1000000000.0) , N_repeat);
    printf("Max Rel Abs Error=%g%%\n\n", MaxAbsRelErr(z__, z4__, NN));

    clock_gettime(CLOCK_MONOTONIC, &start1);     N_iter =5;      funct_cw_newton_iter(Re__, epsD__, z4__, NN, 4, N_repeat);        clock_gettime(CLOCK_MONOTONIC, &finish1);
    printf("Colebrook-White Halley, %d iterations, time= %lf s, rep=%d\n", N_iter, (double)((finish1.tv_sec - start1.tv_sec)+(finish1.tv_nsec - start1.tv_nsec) / 1000000000.0) , N_repeat);
    printf("Max Rel Abs Error=%g%%\n\n", MaxAbsRelErr(z__, z4__, NN));




    //explicit formulas

    //haaland
    clock_gettime(CLOCK_MONOTONIC, &start1);  funct_approx(Re__, epsD__, z_hat, NN, N_repeat,1);  clock_gettime(CLOCK_MONOTONIC, &finish1);
    printf("haaland, time= %lf s, rep=%d\n", (double)((finish1.tv_sec - start1.tv_sec)+(finish1.tv_nsec - start1.tv_nsec) / 1000000000.0)  , N_repeat);
    printf("Max Rel Abs Error=%g%%\n\n", MaxAbsRelErr(z__, z_hat, NN));

    //swamee jain
    clock_gettime(CLOCK_MONOTONIC, &start1);  funct_approx(Re__, epsD__, z_hat, NN, N_repeat, 2);  clock_gettime(CLOCK_MONOTONIC, &finish1);
    printf("swamee jain, time= %lf s, rep=%d\n", (double)((finish1.tv_sec - start1.tv_sec)+(finish1.tv_nsec - start1.tv_nsec) / 1000000000.0)  , N_repeat);
    printf("Max Rel Abs Error=%g%%\n\n", MaxAbsRelErr(z__, z_hat, NN));

    //serghides
    clock_gettime(CLOCK_MONOTONIC, &start1);  funct_approx(Re__, epsD__, z_hat, NN, N_repeat, 3);  clock_gettime(CLOCK_MONOTONIC, &finish1);
    printf("serghides, time= %lf s, rep=%d\n", (double)((finish1.tv_sec - start1.tv_sec)+(finish1.tv_nsec - start1.tv_nsec) / 1000000000.0)  , N_repeat);
    printf("Max Rel Abs Error=%g%%\n\n", MaxAbsRelErr(z__, z_hat, NN));

    //goudar
    clock_gettime(CLOCK_MONOTONIC, &start1);  funct_approx(Re__, epsD__, z_hat, NN, N_repeat, 4);  clock_gettime(CLOCK_MONOTONIC, &finish1);
    printf("goudar-sonnad, time= %lf s, rep=%d\n", (double)((finish1.tv_sec - start1.tv_sec)+(finish1.tv_nsec - start1.tv_nsec) / 1000000000.0) , N_repeat);
    printf("Max Rel Abs Error=%g%%\n\n", MaxAbsRelErr(z__, z_hat, NN));

    //vatan
    clock_gettime(CLOCK_MONOTONIC, &start1);  funct_approx(Re__, epsD__, z_hat, NN, N_repeat, 5);  clock_gettime(CLOCK_MONOTONIC, &finish1);
    printf("vatan 2008, time= %lf s, rep=%d\n", (double)((finish1.tv_sec - start1.tv_sec)+(finish1.tv_nsec - start1.tv_nsec) / 1000000000.0)  , N_repeat);
    printf("Max Rel Abs Error=%g%%\n\n", MaxAbsRelErr(z__, z_hat, NN));

    //praks
    clock_gettime(CLOCK_MONOTONIC, &start1);  funct_approx(Re__, epsD__, z_hat, NN, N_repeat, 6);  clock_gettime(CLOCK_MONOTONIC, &finish1);
    printf("praks     , time= %lf s, rep=%d\n", (double)((finish1.tv_sec - start1.tv_sec)+(finish1.tv_nsec - start1.tv_nsec) / 1000000000.0)  , N_repeat);
    printf("Max Rel Abs Error=%g%%\n\n", MaxAbsRelErr(z__, z_hat, NN));




    //z_prime
    clock_gettime(CLOCK_MONOTONIC, &start1); funct_zp_cw(Re__, z__, zp, NN, 1, N_repeat);    clock_gettime(CLOCK_MONOTONIC, &finish1);
    printf("zp,CW    , time= %lf s, rep=%d\n", (double)((finish1.tv_sec - start1.tv_sec)+(finish1.tv_nsec - start1.tv_nsec) / 1000000000.0)  , N_repeat);



    //if size is 100 then print the results for verification
    
    std::cout << "epsD\tRe\tf \n";
    for (i = 0; i < 10; i++) {
        std::cout << epsD__[i] << "\t" << Re__[i] << "\t" << z__[i] << "\t" << zp[i] <<  "\n";
    }
    std::cout << "\n";


    finish0 = clock();
    printf("\n global timer %f \n\n",  (double)(finish0 - start0)  );

clock_gettime(CLOCK_MONOTONIC, &finish00);
elapsed = (finish00.tv_sec - start00.tv_sec)+(finish00.tv_nsec - start00.tv_nsec) / 1000000000.0;


printf("\n global timer 2= %f",elapsed);

    //destroy dynamic arrays
    delete[] Re__;
    delete[] epsD__;
    delete[] z__;
    delete[] zp;
    delete[] z4__;
    delete[] z_hat;

    std::cout << "Hello World!\n";
}

double err_subs_CW(double* Re__, double* epsD__, double* z_, int NN) {
    int i;
    double temp_ = 0;
    double temp2_ = 0;
    for (i = 0; i < NN; i++) {
        temp2_ = abs(1 / sqrt(z_[i]) + 2 * log10(epsD__[i] / 3.7 + 2.51 / (Re__[i] * sqrt(z_[i]))));
        temp_ = std::max(temp_, temp2_);
    }
    return temp_;
}


double MaxAbsRelErr(double* z_, double* z_hat,int NN) {
    int i;
    double temp_ = 0,temp2_=0;
    for (i = 0; i < NN; i++) {
        temp2_ = abs((z_[i] - z_hat[i]) / z_[i])*100;
        temp_ = std::max(temp_, temp2_);
    }
    return temp_;
}


void ReEpsInitialize(double* Re__, double* epsD__, int N, int type_) {
    int i, j;
    if(type_ == 1){
        for (i = 0; i <= N; i++) {
            for (j = 0; j <= N; j++) {
                epsD__[i + j * (N + 1)] = (1e-6) * pow(0.05 / 1e-6, (double)i / N);  //geometric series
                Re__[i + j * (N + 1)] = (4000) * pow(1e8 / 4000, (double)j / N);  //geometric series
            }
        }
                //111111 1111111111
        epsD__[0]=956423.1564812345561156451;
        Re__[0]=0.0011564564538731584564894;
                  //11111111111111111
    }
    else {
        for (i = 0; i <= N; i++) {
            for (j = 0; j <= N; j++) {
                epsD__[i + j * (N + 1)] = 0 + (0.05-0) *( (double)i / N);  //arithmetic series
                Re__[i + j * (N + 1)] = 4000 + (1e8-4000) * ((double)j / N);  //arithmetic series
            }
        }
        epsD__[0]=956423.1564812345561;
        Re__[0]=0.0011564564538731584;        
    }
}

void funct_cw_iter(double* Re__, double* epsD__, double* z, int NN, int N_iter, int N_repeat) {
    double zz,y1, y2;
    int i, j, k;

    if (N_iter == 3) {

        for (k = 1; k <= N_repeat; k++) {
	#pragma omp parallel
	#pragma omp for
            for (i = 0; i < NN; i++) {
                y1 = epsD__[i] / 3.7;
                y2 = -2.180158299154324 / Re__[i];
                //manual unrolling the loops
                z[i] = 1.325474527619600 * pow(log(y1+y2*log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * (-6)))))), -2); //5iter
            }
        }
    }

    if (N_iter == 4) {

        for (k = 1; k <= N_repeat; k++) {
	#pragma omp parallel
	#pragma omp for
            for (i = 0; i < NN; i++) {
                y1 = epsD__[i] / 3.7;
                y2 = -2.180158299154324 / Re__[i];
                //manual unrolling the loops
                z[i] = 1.325474527619600 * pow(log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * (-6))))), -2); //4iter
            }
        }
    }

    if (N_iter == 5) {

        for (k = 1; k <= N_repeat; k++) {
	#pragma omp parallel
	#pragma omp for
            for (i = 0; i < NN; i++) {
                y1 = epsD__[i] / 3.7;
                y2 = -2.180158299154324 / Re__[i];
                //manual unrolling the loops
                z[i] = 1.325474527619600 * pow(log(y1+y2*log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * (-6)))))), -2); //4iter
            }
        }
    }

    if (N_iter == 30) {

        for (k = 1; k <= N_repeat; k++) {
	#pragma omp parallel
	#pragma omp for
            for (i = 0; i < NN; i++) {
                y1 = epsD__[i] / 3.7;
                y2 = -2.180158299154324 / Re__[i];
                //manual unrolling
		for(j=0;j<8;j++){
	                zz = log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * -6))));//4
		}

                z[i] = 1.325474527619600 / (zz * zz);
            }
        }
    }

    if (N_iter == 50) {
        for (k = 1; k <= N_repeat; k++) {
		#pragma omp parallel
		#pragma omp for
            for (i = 0; i < NN; i++) {
                y1 = epsD__[i] / 3.7;
                y2 = -2.180158299154324 / Re__[i];
                //manual unrolling
		for (j=0;j<13;j++){
	                zz = log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * -6))));//4
		}

                z[i] = 1.325474527619600 / (zz * zz);
            }
        }
    }

}

void funct_zp_cw(double* Re, double* z, double *zp, int NN, int N_iter, int N_repeat) {
    double zz, y1, y2;
    int i, j, k;

    if (N_iter == 1) {
        for (k = 1; k <= N_repeat; k++) {
	#pragma omp parallel
	#pragma omp for
            for (i = 0; i < NN; i++) {
                zp[i] = 1 / (1 / ((-2.51 * 4 / log(10))*z[i] / pow(Re[i],2)*exp((log(10) / 2) / sqrt(z[i]))) - 0.5*Re[i] / z[i]);
            }
        }
    }

}


void funct_cw_newton_iter(double* Re__, double* epsD__, double* z, int NN, int N_iter, int N_repeat) {
    double zz, y1, y2;
    int i, j, k;
    double f_, temp_, fp_, fpp_;
    if (N_iter == 1) {
        for (k = 1; k <= N_repeat; k++) {
	#pragma omp parallel
	#pragma omp for
            for (i = 0; i < NN; i++) {
                //newton-2iter
                y1 = epsD__[i] / 3.7;
                y2 = -2.180158299154324 / Re__[i];
                zz = -6.880288946433447;
                for (j = 0; j < 2; j++) {
                    zz = zz - (zz - log(y1 + y2 * zz)) / (1 - y2 / (y1 + y2 * zz));
                }
                z[i] = 1.325474527619600 * pow(zz, -2);
            }
        }
    }
    if (N_iter == 2) {
        for (k = 1; k <= N_repeat; k++) {
	#pragma omp parallel
	#pragma omp for
            for (i = 0; i < NN; i++) {
                //newton-5iter
                y1 = epsD__[i] / 3.7;
                y2 = -2.180158299154324 / Re__[i];
                zz = -6.880288946433447;
                for (j = 0; j < 5; j++) {
                    zz = zz - (zz - log(y1 + y2 * zz)) / (1 - y2 / (y1 + y2 * zz));
                }
                z[i] = 1.325474527619600 * pow(zz, -2);
            }
        }
    }
    if (N_iter == 3) {
        for (k = 1; k <= N_repeat; k++) {
	#pragma omp parallel
	#pragma omp for
            for (i = 0; i < NN; i++) {
                //halley-2iter
                y1 = epsD__[i] / 3.7;
                y2 = -2.180158299154324 / Re__[i];
                zz = -6.880288946433447;
                for (j = 0; j < 2; j++) {
                    f_ = zz - log(y1 + y2 * zz);
                    temp_ = y2 / (y1 + y2 * zz);
                    fp_ = 1 - temp_;
                    fpp_ = -temp_ * temp_;
                    zz = zz - 2 * f_ * fp_ / (2 * fp_ * fp_ - f_ * fpp_);
                }
                z[i] = 1.325474527619600 * pow(zz, -2);
            }
        }
    }
    if (N_iter == 4) {
        for (k = 1; k <= N_repeat; k++) {
	#pragma omp parallel
	#pragma omp for
            for (i = 0; i < NN; i++) {
                //halley-2iter
                y1 = epsD__[i] / 3.7;
                y2 = -2.180158299154324 / Re__[i];
                zz = -6.880288946433447;
                for (j = 0; j < 5; j++) {
                    f_ = zz - log(y1 + y2 * zz);
                    temp_ = y2 / (y1 + y2 * zz);
                    fp_ = 1 - temp_;
                    fpp_ = -temp_ * temp_;
                    zz = zz - 2 * f_ * fp_ / (2 * fp_ * fp_ - f_ * fpp_);
                }
                z[i] = 1.325474527619600 * pow(zz, -2);
            }
        }
    }

}


void funct_approx(double* Re, double* epsD, double* z, int NN, int N_repeat, int method_) {
    int i, j, k;
    double A, B, C,y,y1,y2; //sergh
    double aa, bb, dd, ss, qq, gg, zz, dla, dcfa; //goudar
    
    if (method_ == 1) { //haland 
        for (k = 1; k <= N_repeat; k++) {
	#pragma omp parallel
	#pragma omp for
            for (i = 0; i < NN; i++) {
                z[i] = pow(-1.8 * log10(pow(epsD[i] / 3.7, 1.11) + 6.9 / Re[i]),-2);
            }
        }
    }

    if (method_ == 2) { //swammee-jain 
        for (k = 1; k <= N_repeat; k++) {
	#pragma omp parallel
	#pragma omp for
            for (i = 0; i < NN; i++) {
                z[i] = 0.25*pow(log10(epsD[i] / 3.7 + 5.074 / pow(Re[i],0.9)),-2);
            }
        }
    }

    if (method_ == 3) { //serghides
        for (k = 1; k <= N_repeat; k++) {
	#pragma omp parallel
	#pragma omp for
            for (i = 0; i < NN; i++) {
                y1 = epsD[i] / 3.7;
                y2 = 2.51 / Re[i];
                A = -0.868588963806504 * log(y1 + 4.7808764940239043824701195219124*y2);
                B = -0.868588963806504 * log(y1 + y2*A );
                C = -0.868588963806504 * log(y1 + y2 * B );
                z[i] =pow( A - (B - A) * (B - A) / (C - 2 * B + A),-2);
            }
        }
    }

    if (method_ == 4) { //goudar-sonnad
        for (k = 1; k <= N_repeat; k++) {
	#pragma omp parallel
	#pragma omp for
            for (i = 0; i < NN; i++) {
                aa = 2 / log(10);
                bb = epsD[i] / 3.7;
                dd = log(10) * Re[i] / 5.02;
                ss = bb * dd + log(dd);
                qq = pow(ss, ss / (ss + 1));
                gg = bb * dd + log(dd / qq);
                zz = log(qq / gg);
                dla = zz * gg / (gg + 1);
                dcfa = dla * (1 + zz / 2 / ((gg + 1) * (gg + 1) + zz / 3 * (2 * gg - 1)));
                z[i] = pow(aa * (log(dd / qq) + dcfa),-2);
                //z[i] = 1 / (z[i] * z[i]);
            }
        }
    }

    if (method_ == 5) { //vatan 2008
        for (k = 1; k <= N_repeat; k++) {
	#pragma omp parallel
	#pragma omp for
            for (i = 0; i < NN; i++) {
                ss = 0.124 * Re[i] * epsD[i] + log(0.4587 * Re[i]);
                z[i] = pow(0.8686 * log(0.4587 * Re[i] / pow(ss - 0.31, ss / (ss + 0.9633))), -2);
                //z[i] = 1 / (z[i] * z[i]);
            }
        }
    }


    if (method_ == 6) { //praks
        for (k = 1; k <= N_repeat; k++) {
	#pragma omp parallel
	#pragma omp for
            for (i = 0; i < NN; i++) {
                ss = 74205.5 + 1000. * epsD[i]*Re[i];
                A = 8 + (2 / log(10))*log(16. / Re[i] + epsD[i] / 3.7);
                B = -74914381.46*pow(ss,-2);
                C = 1391459721232.67 * pow(ss, -3);
                zz = 8 - A - 0.5*A*A*B;
                zz = (-2 / log(10))*log(2.51 / Re[i] *zz + epsD[i] / 3.7);
                zz = (-2 / log(10))*log(2.51 / Re[i] *zz + epsD[i] / 3.7);
                zz = (-2 / log(10))*log(2.51 / Re[i] *zz + epsD[i] / 3.7);
                z[i] = pow(zz, -2);
            }
        }
    }
    

}





//intel i7 3610qm windows 10 cygwin
//matrix size=3164*3164 or 10010896 elements
//initiliaze Re & epsD duration= 0.503083 s
//Colebrook-White, 50 iterations, time= 1.376782 s, rep=10
//err substituting in CW equation = 0.000540684
//Colebrook-White, 3 iterations, time= 1.711198 s, rep=10
//Max Rel Abs Error=0.0157835%
//Colebrook-White, 4 iterations, time= 1.372159 s, rep=10
//Max Rel Abs Error=0%
//Colebrook-White, 5 iterations, time= 1.676324 s, rep=10
//Max Rel Abs Error=0.0157835%
//Colebrook-White, 30 iterations, time= 1.336063 s, rep=10
//Max Rel Abs Error=0%
//Colebrook-White Newton, 2 iterations, time= 0.933774 s, rep=10
//Max Rel Abs Error=0.0139539%
//Colebrook-White Newton, 5 iterations, time= 2.261522 s, rep=10
//Max Rel Abs Error=0.0139595%
//Colebrook-White Hallay, 2 iterations, time= 1.003900 s, rep=10
//Max Rel Abs Error=0.0139165%
//Colebrook-White Halley, 5 iterations, time= 2.402742 s, rep=10
//Max Rel Abs Error=0.0139595%
//haaland, time= 1.140690 s, rep=10
//Max Rel Abs Error=1.42789%
//swamee jain, time= 1.164965 s, rep=10
//Max Rel Abs Error=3.66134%
//serghides, time= 1.106395 s, rep=10
//Max Rel Abs Error=0.0167052%
//goudar-sonnad, time= 1.913201 s, rep=10
//Max Rel Abs Error=0.0139595%
//vatan 2008, time= 1.317549 s, rep=10
//Max Rel Abs Error=0.0539776%
//praks     , time= 1.502662 s, rep=10
//Max Rel Abs Error=0.0222661%
//zp,CW    , time= 0.716674 s, rep=10
// global timer 177249.000000
// global timer 2= 24.154802Hello World!
