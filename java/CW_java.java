// import java.util.Scanner;
import java.io.*;


import java.util.concurrent.*;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Random;
import java.lang.*;

class CW_java {

    public static void main (String[] args) {

        long startTime;
        long estimatedTime;
        int n;

        // create a scanner so we can read the command-line input
        //Scanner scanner = new Scanner(System.in);

        // create instance of Random class
        Random rand = new Random();

        //n = Integer.parseInt(args[0]);
        n=9000;
        n=10;
        n=50;
        n=10;
        n=20;
        n=10000000;
        //n=1000000;
        
        int i,j,k;

        
        
        double Re[]=new double[n];
        double epsD[]=new double[n];
        double z_[]=new double[n];
        double z_hat[]=new double[n];
        double z,y1,y2,y3;
        
        
        for (i=0;i<n;i++) {
                Re[i] = 4000*Math.pow(1e8/4000,rand.nextDouble());
                epsD[i]= (1e-6)*Math.pow(0.05/1e-6,rand.nextDouble());
        }

        System.out.println("Copy Matrix: "+n);

        startTime = System.currentTimeMillis();


        for(i=0;i<n;i++) {
        	y1 = epsD[i] / 3.7;
        	y2 = -2.180158299154324 / Re[i];
        	z=-8;
        	for(j=0;j<30;j++) {
        		z=Math.log(y1+y2*z);
        	}
        	//z=Math.log(y1+y2*Math.log(y1 + y2 * Math.log(y1 + y2 * Math.log(y1 + y2 * Math.log(y1 + y2 * (-6))))));
        	z_[i] = 1.325474527619600 /(z*z);
        }
        
        
        estimatedTime = System.currentTimeMillis() - startTime;
        NumberFormat formatter = new DecimalFormat("#00000.00000000");
        System.out.println("    Time CW 5 iter: " + formatter.format(estimatedTime  / 1000d) + " seconds");
        
        
        
        if(n<20) {
        	for(i=0;i<n;i++) {
        		System.out.println(Re[i]+","+epsD[i]+","+z_[i]);
        	}
        }

        
        
    }

}



