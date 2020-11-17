
//intelc++, windows 10, core i7 3610 qm


#include <iostream>
#include <chrono>
#include <sys/types.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>


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

    N_repeat = 10;
 

    std::cout.precision(17);

    clock_t start, finish;
    srand(time(NULL));

    printf("matrix size=%d*%d or %d elements\n\n", N + 1, N + 1, NN);

    //initialize Re & epsD
    start = clock();
    ReEpsInitialize(Re__, epsD__, N, geometric_series); 
    finish = clock();
    printf("initiliaze Re & epsD duration= %lf s\n\n",  (double)(finish - start) / CLOCKS_PER_SEC );


    //cw 50 iter
    start = clock();       N_iter = 50;       funct_cw_iter(Re__, epsD__, z__, NN, N_iter,N_repeat);        finish = clock();
    printf("Colebrook-White, %d iterations, time= %lf s, rep=%d\n", N_iter, (double)(finish - start) / CLOCKS_PER_SEC,N_repeat);
    printf("err substituting in CW equation = %g\n\n", err_subs_CW(Re__, epsD__, z__, NN));




    //cw 3,4,5,30 iters
    start = clock();       N_iter = 3;      funct_cw_iter(Re__, epsD__, z4__, NN, N_iter, N_repeat);       finish = clock();
    printf("Colebrook-White, %d iterations, time= %lf s, rep=%d\n", N_iter, (double)(finish - start) / CLOCKS_PER_SEC , N_repeat);
    printf("Max Rel Abs Error=%g%%\n\n", MaxAbsRelErr(z__, z4__, NN));

    start = clock();       N_iter = 4;        funct_cw_iter(Re__, epsD__, z4__, NN, N_iter, N_repeat);        finish = clock();
    printf("Colebrook-White, %d iterations, time= %lf s, rep=%d\n", N_iter, (double)(finish - start) / CLOCKS_PER_SEC, N_repeat);
    printf("Max Rel Abs Error=%g%%\n\n", MaxAbsRelErr(z__, z4__, NN));

    start = clock();       N_iter = 5;        funct_cw_iter(Re__, epsD__, z4__, NN, N_iter, N_repeat);       finish = clock();
    printf("Colebrook-White, %d iterations, time= %lf s, rep=%d\n", N_iter, (double)(finish - start) / CLOCKS_PER_SEC, N_repeat);
    printf("Max Rel Abs Error=%g%%\n\n", MaxAbsRelErr(z__, z4__, NN));

    start = clock();      N_iter = 30;        funct_cw_iter(Re__, epsD__, z4__, NN, N_iter, N_repeat);       finish = clock();
    printf("Colebrook-White, %d iterations, time= %lf s, rep=%d\n", N_iter, (double)(finish - start) / CLOCKS_PER_SEC, N_repeat);
    printf("Max Rel Abs Error=%g%%\n\n", MaxAbsRelErr(z__, z4__, NN));



    //newton & halley

    start = clock();        N_iter = 2;      funct_cw_newton_iter(Re__, epsD__, z4__, NN, 1, N_repeat);        finish = clock();
    printf("Colebrook-White Newton, %d iterations, time= %lf s, rep=%d\n", N_iter, (double)(finish - start) / CLOCKS_PER_SEC , N_repeat);
    printf("Max Rel Abs Error=%g%%\n\n", MaxAbsRelErr(z__, z4__, NN));

    start = clock();      N_iter = 5;       funct_cw_newton_iter(Re__, epsD__, z4__, NN, 2, N_repeat);       finish = clock();
    printf("Colebrook-White Newton, %d iterations, time= %lf s, rep=%d\n", N_iter, (double)(finish - start) / CLOCKS_PER_SEC, N_repeat);
    printf("Max Rel Abs Error=%g%%\n\n", MaxAbsRelErr(z__, z4__, NN));

    start = clock();    N_iter = 2;      funct_cw_newton_iter(Re__, epsD__, z4__, NN, 3, N_repeat);        finish = clock();
    printf("Colebrook-White Hallay, %d iterations, time= %lf s, rep=%d\n", N_iter, (double)(finish - start) / CLOCKS_PER_SEC, N_repeat);
    printf("Max Rel Abs Error=%g%%\n\n", MaxAbsRelErr(z__, z4__, NN));

    start = clock();     N_iter =5;      funct_cw_newton_iter(Re__, epsD__, z4__, NN, 4, N_repeat);        finish = clock();
    printf("Colebrook-White Halley, %d iterations, time= %lf s, rep=%d\n", N_iter, (double)(finish - start) / CLOCKS_PER_SEC, N_repeat);
    printf("Max Rel Abs Error=%g%%\n\n", MaxAbsRelErr(z__, z4__, NN));




    //explicit formulas

    //haaland
    start = clock();  funct_approx(Re__, epsD__, z_hat, NN, N_repeat,1);  finish = clock();
    printf("haaland, time= %lf s, rep=%d\n", (double)(finish - start) / CLOCKS_PER_SEC , N_repeat);
    printf("Max Rel Abs Error=%g%%\n\n", MaxAbsRelErr(z__, z_hat, NN));

    //swamee jain
    start = clock();  funct_approx(Re__, epsD__, z_hat, NN, N_repeat, 2);  finish = clock();
    printf("swamee jain, time= %lf s, rep=%d\n", (double)(finish - start) / CLOCKS_PER_SEC , N_repeat);
    printf("Max Rel Abs Error=%g%%\n\n", MaxAbsRelErr(z__, z_hat, NN));

    //serghides
    start = clock();  funct_approx(Re__, epsD__, z_hat, NN, N_repeat, 3);  finish = clock();
    printf("serghides, time= %lf s, rep=%d\n", (double)(finish - start) / CLOCKS_PER_SEC , N_repeat);
    printf("Max Rel Abs Error=%g%%\n\n", MaxAbsRelErr(z__, z_hat, NN));

    //goudar
    start = clock();  funct_approx(Re__, epsD__, z_hat, NN, N_repeat, 4);  finish = clock();
    printf("goudar-sonnad, time= %lf s, rep=%d\n", (double)(finish - start) / CLOCKS_PER_SEC, N_repeat);
    printf("Max Rel Abs Error=%g%%\n\n", MaxAbsRelErr(z__, z_hat, NN));

    //vatan
    start = clock();  funct_approx(Re__, epsD__, z_hat, NN, N_repeat, 5);  finish = clock();
    printf("vatan 2008, time= %lf s, rep=%d\n", (double)(finish - start) / CLOCKS_PER_SEC , N_repeat);
    printf("Max Rel Abs Error=%g%%\n\n", MaxAbsRelErr(z__, z_hat, NN));

    //praks
    start = clock();  funct_approx(Re__, epsD__, z_hat, NN, N_repeat, 6);  finish = clock();
    printf("praks     , time= %lf s, rep=%d\n", (double)(finish - start) / CLOCKS_PER_SEC , N_repeat);
    printf("Max Rel Abs Error=%g%%\n\n", MaxAbsRelErr(z__, z_hat, NN));




    //z_prime
    start = clock(); funct_zp_cw(Re__, z__, zp, NN, 1, N_repeat);    finish = clock();
    printf("zp,CW    , time= %lf s, rep=%d\n", (double)(finish - start) / CLOCKS_PER_SEC , N_repeat);



    //if size is 100 then print the results for verification
    if (N <= 10) {
        std::cout << "epsD\tRe\tf \n";
        for (i = 0; i < NN; i++) {
            std::cout << epsD__[i] << "\t" << Re__[i] << "\t" << z__[i] << "\t" << zp[i] <<  "\n";
        }
        std::cout << "\n";
    }



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
    }
    else {
        for (i = 0; i <= N; i++) {
            for (j = 0; j <= N; j++) {
                epsD__[i + j * (N + 1)] = 0 + (0.05-0) *( (double)i / N);  //arithmetic series
                Re__[i + j * (N + 1)] = 4000 + (1e8-4000) * ((double)j / N);  //arithmetic series
            }
        }
    }
}

void funct_cw_iter(double* Re__, double* epsD__, double* z, int NN, int N_iter, int N_repeat) {
    double zz,y1, y2;
    int i, j, k;

    if (N_iter == 3) {
        for (k = 1; k <= N_repeat; k++) {
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
            for (i = 0; i < NN; i++) {
                y1 = epsD__[i] / 3.7;
                y2 = -2.180158299154324 / Re__[i];
                //manual unrolling
                zz = log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * -6))));//4
                zz = log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * zz))));//8
                zz = log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * zz))));//12
                zz = log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * zz))));//16
                zz = log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * zz))));//18
                zz = log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * zz))));//24
                zz = log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * zz))));//28
                zz = log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * zz))));//32

                z[i] = 1.325474527619600 / (zz * zz);
            }
        }
    }

    if (N_iter == 50) {
        for (k = 1; k <= N_repeat; k++) {
            for (i = 0; i < NN; i++) {
                y1 = epsD__[i] / 3.7;
                y2 = -2.180158299154324 / Re__[i];
                //manual unrolling
                zz = log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * -6))));//4
                zz = log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * zz))));//8
                zz = log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * zz))));//12
                zz = log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * zz))));//16
                zz = log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * zz))));//18
                zz = log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * zz))));//24
                zz = log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * zz))));//28
                zz = log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * zz))));//32
                zz = log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * zz))));//36
                zz = log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * zz))));//40
                zz = log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * zz))));//44
                zz = log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * zz))));//48
                zz = log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * log(y1 + y2 * zz))));//52

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
            for (i = 0; i < NN; i++) {
                z[i] = pow(-1.8 * log10(pow(epsD[i] / 3.7, 1.11) + 6.9 / Re[i]),-2);
            }
        }
    }

    if (method_ == 2) { //swammee-jain 
        for (k = 1; k <= N_repeat; k++) {
            for (i = 0; i < NN; i++) {
                z[i] = 0.25*pow(log10(epsD[i] / 3.7 + 5.074 / pow(Re[i],0.9)),-2);
            }
        }
    }

    if (method_ == 3) { //serghides
        for (k = 1; k <= N_repeat; k++) {
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
            for (i = 0; i < NN; i++) {
                ss = 0.124 * Re[i] * epsD[i] + log(0.4587 * Re[i]);
                z[i] = pow(0.8686 * log(0.4587 * Re[i] / pow(ss - 0.31, ss / (ss + 0.9633))), -2);
                //z[i] = 1 / (z[i] * z[i]);
            }
        }
    }


    if (method_ == 6) { //praks
        for (k = 1; k <= N_repeat; k++) {
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

