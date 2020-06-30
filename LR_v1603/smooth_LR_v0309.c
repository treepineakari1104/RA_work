// 筋活動を 1. (整流 + Low-pass) でスムージング. 
// コンパイル: gcc -lfftw3 smooth_LR.c -Wall
// 印刷: OK (15/6/23)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>  // FFTW 

#define    MAX        100   //  ピーク数
#define     DT     0.0001   //  サンプリング [s]


// 信号 x[t] (サイズ: num) を平滑化.
void smooth_xt(int n, double f_m, double x[], double x_sm[] );
// Detect peak in smoothed motor activity. 
int calc_peak(int n, double x[], int peak[]);

int main( int argc, char *argv[] )
{    
  if( argc != 4 )   
    {  printf("Usage:./a.out *f_in $num *f_out\n");  exit(1);  }
  
  FILE *f_in, *f_out, *f_peak; 

  if( ( f_in = fopen(argv[1], "r"))==NULL )  // Input: Voltage, L, R
    {  printf("file \"%s\" is not found\n",argv[1]);   exit(1);  }
  if( ( f_out = fopen(argv[3], "w"))==NULL ) // Output: Smoothed L R
    {  printf("file \"%s\" is not found\n",argv[3]);   exit(1);  }
  
  int    i, num= atoi(argv[2]);    
  double tmp, *L, *R, *L_sm, *R_sm, ave_L= 0, ave_R= 0;
  // Assign  Memory
  L= (double*)calloc( num, sizeof(double) );
  R= (double*)calloc( num, sizeof(double) );
  L_sm= (double*)calloc( num, sizeof(double) );
  R_sm= (double*)calloc( num, sizeof(double) );
  
  for ( i=0; i<num; i++ ) // Read data: Volt R L
    {  fscanf(f_in, "%lf %lf %lf", &tmp, &R[i], &L[i]);    }
  
  double f_lp= 2.0;    // Low-pass Filter の振動数.
  smooth_xt( num, f_lp, L, L_sm);  // Low-pass Filter: L
  smooth_xt( num, f_lp, R, R_sm);  // Low-pass Filter: R
  
  for (i=0; i<num; i++)
  {   ave_L += L[i];   ave_R += R[i];   }
  ave_L= ave_L/num;     ave_R= ave_R/num;
    
  for (i=0; i<num; i++)   
    {  
      if ( i%10==0 )
	{  fprintf( f_out, "%2.3f  %f %f %f %f\n", i* 0.0001, R[i]-ave_R,
               L[i]-ave_L, R_sm[i], L_sm[i] );  }
    } 
  
  f_peak = fopen("peak.dat","w");
  if ( f_peak == NULL)  {  printf("Cannot open *f_peak\n");  exit(1);  }
  
  int  n_p, peak[MAX]= {0};
  n_p= calc_peak( num, R_sm, peak);  
  fprintf( f_peak, "%d \n", n_p);   
  for (i=0; i<n_p; i++)
    {  
      peak[i]=  (int)( (peak[i]+5)/10)* 10;  
      fprintf( f_peak, "%d %f \n", peak[i], R_sm[ peak[i] ]);         
    }   
  
  fclose(f_in);     fclose(f_out);     fclose(f_peak);    
  free(L);    free(R);    free(L_sm);    free(R_sm);
  return 0;
}


// Smoothing the muscle potential (L: Left, R: Right)
void smooth_xt(int n, double f_m, double x[], double x_sm[] )
{
  fftw_complex	*in, *out;              fftw_plan  plan;
  int  i, n_lp= (int)(f_m*n*DT);        double    ave= 0; 
  
  // メモリの確保
  in   = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);
  out  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);  
  plan = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);      
  
  for (i=0; i<n; i++)  // 信号の読み込み
    {  in[i][0]= x[i];   in[i][1] = 0.0;   ave += in[i][0];  }
  
  ave= ave/n;          // 前処理: 全波整流
  for (i=0; i<n; i++)  {  in[i][0]= fabs( in[i][0]-ave );    }    
  
  ave= 0;              // 前処理: 全波整流後の平均値の計算. 
  for (i=0; i<n; i++)  {  ave += in[i][0];    }
  
  ave= ave/n;          // 前処理: 信号の平均値を引く
  for (i=0; i<n; i++)  {  in[i][0] -= ave;    }
  
  fftw_execute(plan);      // フーリエ変換  
  fftw_destroy_plan(plan);
    
  for (i=0; i<n; i++)  // 高周波数のスペクトルを 0 にする.
    { 
      in[i][0]= out[i][0]/n;    in[i][1]= out[i][1]/n; 
      //  x_j, x_{n-j}: 同じ周波数であることに注意！
      if ( n_lp< i &&  i< n-n_lp )  {  in[i][0]= 0.0;   in[i][1]= 0.0;  }
    }
  in[0][0]= 0;    in[0][1]= 0;  // 定数部分: 0 にする。
    
  plan = fftw_plan_dft_1d(n, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(plan);     // フーリエ逆変換    
  for (i= 0; i<n; i++)    // スムージングデータの格納
    {  x_sm[i]= out[i][0];    }
  
  fftw_free(in);   fftw_free(out);     fftw_destroy_plan(plan);
}


int calc_peak(int n, double x[], int peak[])
{ // 平滑化された信号のピーク時刻を検出
  int i, tmp= 0, info= 1;
  
  for (i=500; i<n-2; i++)
    {
      if ( x[i]<0 && 0<= x[i+1])  // ゼロ交差判定
	{  info= 1;  }	  
      
      if ( x[i]<= x[i+1] && x[i+1] >= x[i+2] && info==1 )  // ピーク検出
	{  peak[tmp]= i+1;  info= 0;  tmp++;  }
    }
  return tmp;  // ピーク数を返す.
}
