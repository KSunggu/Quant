#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "Koba_six.h"

int it, nx, ny, m;
double h, h2, eps, pi, cool, hot, D, K, alpha, beta, lam, dt, xleft, xright, 
        yleft, yright, e4, **work, seed_r;

int main()
{  
    extern int it, nx, ny, m;
    extern double h, h2, eps, pi, cool, hot, D, K, alpha, beta, lam, dt, 
            xleft, xright, yleft, yright, e4, **work, seed_r;    
    int ii, jj, x1, x2, y1, y2, i, j, max_it, ns, count = 1;
    double t, deno, phix, phiy, gphi, gphi2, gphi6, T,
            **flag, **oT, **nT, **oc, **nc, **W, **W2, **interx, **intery, 
            **nc1;
    double elapsed;
    clock_t start, end;
    FILE *fp, *my;
    
    start = clock();
    
    pi = 4.0*atan(1.0);
    
    nx = gnx;  ny = gny;

    h = 0.4;   h2 = h*h;

    xleft = 0.0, xright = h*nx;
    yleft = 0.0, yright = h*ny;

    seed_r = 2.0*3.462;
    
    ////////////////////////////////////////////////// Sample ///////////////////////////////////////////////////////

    nc1 = dmatrix(0, nx, 0, ny);
    oT = dmatrix(0, nx, 0, ny);
    nT = dmatrix(0, nx, 0, ny);
    oc = dmatrix(0, nx, 0, ny);
    nc = dmatrix(0, nx, 0, ny);
    work = dmatrix(0, nx, 0, ny);
    W = dmatrix(0, nx, 0, ny);
    W2 = dmatrix(0, nx, 0, ny);
    interx = dmatrix(0, nx, 0, ny);
    intery = dmatrix(0, nx, 0, ny);
    flag = dmatrix(0, nx, 0, ny);

    initial(oT, oc);
    mat_copy(nT, oT, 0, nx, 0, ny);
    mat_copy(nc, oc, 0, nx, 0, ny);
    mat_copy(nc1, oc, 0, nx, 0, ny);
    print_data(nT, nc, count);

 
    // flag index
//     x1 = (int)((xright*0.5 - seed_r)/h);     x2 = nx-x1;
//     y1 = (int)((yright*0.5 - seed_r)/h);     y2 = ny-y1;
       x1 = 1; x2 = nx - 1; 
       y1 = 1; y2 = ny - 1; 
    t = 0.0;   
     for (it=1; it<=max_it; it++){
        
        if (it%1000 == 0) 
            printf("iteration= %d\n", it);

        t += dt;
        zero_matrix(flag, x1, x2, y1, y2);
 
        
        //* flag  setting *//
 ////////////////////////////////////////////////////// Sample ///////////////////////////////////////////////////////////

  
	   for (i=1; i<=nx/2; i++) {        
           if (oc[nx/2-i+1][ny/2] < 0.01) {
               x1 = nx/2-i+1-(2*m+1); // 0.01보다 큰 쪽으로 5포인트 앞
               i = nx; 
           } 
       }
       x2 = nx - x1 + 1;
       y1 = x1;   y2 = x2;
       
       newiter {
           if ( flag[i][j] > 0.5 ) {
               phix = (oc[i+1][j] - oc[i-1][j])/(2.0*h);
               phiy = (oc[i][j+1] - oc[i][j-1])/(2.0*h);
               gphi = phix*phix + phiy*phiy; // |\nabla phi|^2
               gphi2 = pow(gphi, 2) + eps;
               gphi6 = pow(gphi, 3) + eps;
               
 ////////////////////////////////////////////////////// Sample ///////////////////////////////////////////////////////////

      

           }           
       }
       
    
       newiter {
           if (flag[i][j]*flag[i-1][j]*flag[i+1][j]*flag[i][j-1]*flag[i][j+1]>0.5) {             
               //  step 1
 ////////////////////////////////////////////////////// Sample ///////////////////////////////////////////////////////////
               //  step 2
 ////////////////////////////////////////////////////// Sample ///////////////////////////////////////////////////////////

               
           }
       }
       
       // step 3
       ijloopin {
           ////////////////////////////////////////////////////// Sample ///////////////////////////////////////////////////////////
       }
       
       mat_copy(oT, nT, 1, nx-1, 1, ny-1);  
       mat_copy(oc, nc, 1, nx-1, 1, ny-1);  
      
    
 
      
//         printf("iteration= %d\n", it);
//         if (it==max_it){
//         my = fopen("test3.m","w");
//         fprintf(my,"T=[ ");
//         for (j=0; j<=ny; j++)
//         { for (i=0; i<=nx; i++)
//         { fprintf(my," %f ",flag[i][j]); }
//         fprintf(my,"\n");}
//         fprintf(my,"]; mesh(T')");
//         fclose(my); }
           
      if( (it % ns)==0) {   
	        count++;	 
            print_data(nT,nc,count);
            printf("print out counts %d\n",count); 
            end = clock();
            elapsed = ((float) (end - start)) / CLOCKS_PER_SEC;
            my = fopen("data/remarks.m","a");
            fprintf(my, "\n  counts=    %d;     CPU_Time=   %f;\n", count-1, elapsed);
            fclose(my); 
            
//                my = fopen("test.m","w");
//   fprintf(my,"close all; T=[ ");
//   for (j=1; j<=ny; j++)
//   { for (i=1; i<=nx; i++)
//       { fprintf(my," %f ",oT[i][j]); }
//   fprintf(my,"\n");}
//   fprintf(my,"]; mesh(T');");
//   fclose(my); 
//   exit(0);
            
         } 
	   
       
     } //end-iteration 
  
   fp = fopen("cry.m","a");
   fclose(fp);
   end = clock();
   elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
   printf("Time elapsed %f\n",elapsed);
   my = fopen("data/remarks.m","a");
   fprintf(my, "\n Finally_time_elapsed = %f;\n", elapsed);
   fclose(my);

   return 0;
}


void initial(double **T, double **c)
{
   extern int nx, ny;
   extern double cool, hot, xright, yright, seed_r ;
   int i, j;
   double rp, x, y;
   FILE *my;
   
   ijloop {

      x = ((double)i) * h;
      y = ((double)j) * h;

      rp = sqrt( pow(x-xright/2.0,2) + pow(y-yright/2.0,2) );
      c[i][j] = 0.5*(1.0 + tanh((seed_r - rp)/(1.0*sqrt(2))));
      T[i][j] = cool;
      if (c[i][j] > 0.5)
          T[i][j] = hot;
    }
//    my = fopen("test.m","w");
//   fprintf(my,"close all; T=[ ");
//   for (j=1; j<=ny; j++)
//   { for (i=1; i<=nx; i++)
//       { fprintf(my," %f ",c[i][j]); }
//   fprintf(my,"\n");}
//   fprintf(my,"]; mesh(T'); axis image");
//   fclose(my); 
//   exit(0);
   
   
}

void zero_matrix(double **a, int xl, int xr, int yl, int yr)
{
   int i, j;

   for (i=xl; i<=xr; i++)
      for (j=yl; j<=yr; j++){
        a[i][j] = 0.0;	 
  } 
}

void mat_copy(double **a, double **b, int xl, int xr, int yl, int yr)
{
    int i, j;
    
    for (i=xl; i<=xr; i++)
        for (j=yl; j<=yr; j++)
            a[i][j] = b[i][j];
}

void print_mat(FILE *fptr, double **a, int nrl, int nrh, int ncl, int nch)
{
    int i, j;
    
    for(j = ncl; j <= nch; j++) {
        for(i = nrl; i <= nrh; i++)
            fprintf(fptr, " %f", a[i][j]);
        
        fprintf(fptr, "\n");
    }

}


void print_data(double **t, double **c, int count)
{
    extern int nx, ny;
    char buffer[20], buffert[20];
    FILE *ft, *fc;
    
    sprintf(buffer,"data/cry%d.m",count);
    sprintf(buffert,"data/temp%d.m",count);
    fc = fopen(buffer,"w");
    fclose(fc);
    fc = fopen(buffer,"a"); 
    ft = fopen(buffert,"w");
    fclose(ft);
    ft = fopen(buffert,"a");
    
    print_mat(ft, t, 0, nx, 0, ny);
    print_mat(fc, c, 0, nx, 0, ny);
    
    fclose(ft);
    fclose(fc);
}

double F(double c)
{
	double value;

	//value = 0.25*pow(c*c-1.0,2);
    //value = pow(c*(1-c),1); 
    value = pow(c*(c-1),2);
	return value;
}


double *dvector (long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
   double *v;

   v=(double *) malloc((nh-nl+1+1)*sizeof(double));
   return v-nl+1;

}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
{
   double **m;
   long i, nrow=nrh-nrl+1+1, ncol=nch-ncl+1+1;

   m=(double **) malloc((nrow)*sizeof(double*));
   m+=1;
   m-=nrl;

   m[nrl]=(double *) malloc((nrow*ncol)*sizeof(double));
   m[nrl]+=1;

   m[nrl]-=ncl;

   for (i=nrl+1; i<=nrh; i++) m[i]=m[i-1]+ncol;

   return m;
}
