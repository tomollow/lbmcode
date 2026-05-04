// Poiseuille flow 
// Lattice Boltzmann Method.
// The source code is written in C programming language.
// Copyright @2021 Takeshi Seta All Rights Reserved.
//
// 6  2  5
//    |
// 3--0--1
//    |
// 7  4  8
//
// rho: density
// u  : horizontal velocity component (time = n)
// v  : vertical velocity component (time = n)
// un : horizontal velocity component (time = n -1)
// vn : vertical velocity component (time = n - 1)
// gx : body force (x-axis)
// gy : body force (y-axis)
// h  : channel width
// nx : number of grid points (x-axis)
// ny : number of grid points (y-axis)
// f  : distribution function
// f0 : equilibrium distribution function
// cx : discrete velocity (x-axis)
// cy : discrete velocity (y-axis)
// tau: relaxation time
// nu : kinematic viscosity

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

# define DIM 64

int main(int argc, char *argv[])
{
  FILE   *fp;
  int    nx = 20, ny = 20, time = 0, loop1, loop2;
  int    i, j, k, in, jn, nb;
  double u[DIM][DIM], v[DIM][DIM], un[DIM][DIM], vn[DIM][DIM], rho[DIM][DIM];
  double f[9][DIM][DIM], f0[9][DIM][DIM], ftmp[9][DIM][DIM], cx[9], cy[9];
  double gx = 0.00001, gy = 0.0, tmp, u2, h, nu, norm, tau = 0.56;
  char   a[DIM][DIM];

  // tau はコマンドライン引数で上書きできるようにする。
  if(argc >= 2){
    tau = atof(argv[1]);
  }
  nu = (tau - 0.5)/3.0;
  h = (double)ny;

  // 密度 1、速度 0 の静止場から計算を始める。
  for(i = 0; i <= nx; i++){ for(j = 0; j <= ny; j++){
    u[i][j] = 0.0; v[i][j] = 0.0; un[i][j] = 0.0; vn[i][j] = 0.0; rho[i][j] = 1.0;
  } }

  cx[0] = 0.0; cy[0] = 0.0;
  cx[1] = 1.0; cy[1] = 0.0;  cx[2] = 0.0; cy[2] = 1.0;
  cx[3] =-1.0; cy[3] = 0.0;  cx[4] = 0.0; cy[4] =-1.0;
  cx[5] = 1.0; cy[5] = 1.0;  cx[6] =-1.0; cy[6] = 1.0;
  cx[7] =-1.0; cy[7] =-1.0;  cx[8] = 1.0; cy[8] =-1.0;

  for(i = 0; i <= nx; i++){ for(j = 0; j <= ny; j++){
    u2 = u[i][j]*u[i][j] + v[i][j]*v[i][j];      
    f0[0][i][j] = rho[i][j]*(1.0 -3.0/2.0*u2)*4.0/9.0;
    for(k = 1; k <= 4; k++){
      tmp = cx[k]*u[i][j] + cy[k]*v[i][j];      
      f0[k][i][j] = rho[i][j]*(1.0 +3.0*tmp +9.0/2.0*tmp*tmp -3.0/2.0*u2)/9.0;
    }
    for(k = 5; k <= 8; k++){
      tmp = cx[k]*u[i][j] + cy[k]*v[i][j];      
      f0[k][i][j] = rho[i][j]*(1.0 +3.0*tmp +9.0/2.0*tmp*tmp -3.0/2.0*u2)/36.0;
    }
  } }

  for(i = 0; i <= nx; i++){ for(j = 0; j <= ny; j++){ for(k = 0; k <= 8; k++){
    f[k][i][j] = f0[k][i][j];
  } } }

  // 定常解に近づくまで、LBM を反復して速度場を更新する。
  for(loop1 = 0; loop1 < 500; loop1++){
  for(loop2 = 0; loop2 < 200; loop2++){
    time++;

    // 収束判定のため、1 ステップ前の速度場を保存する。
    for(i = 0; i <= nx; i++){ for(j = 0; j <= ny; j++){
      un[i][j] = u[i][j]; vn[i][j] = v[i][j];
    } }

    // BGK 衝突モデルで分布関数を平衡分布へ緩和させる。
    for(i = 0; i <= nx; i++){ for(j = 0; j <= ny; j++){
      u2 = u[i][j]*u[i][j] + v[i][j]*v[i][j];      
      f0[0][i][j] = rho[i][j]*(1.0 -3.0/2.0*u2)*4.0/9.0;
      for(k = 1; k <= 4; k++){
        tmp = cx[k]*u[i][j] + cy[k]*v[i][j];      
        f0[k][i][j] = rho[i][j]*(1.0 +3.0*tmp +9.0/2.0*tmp*tmp -3.0/2.0*u2)/9.0;
      }  
      for(k = 5; k <= 8; k++){
        tmp = cx[k]*u[i][j] + cy[k]*v[i][j];      
        f0[k][i][j] = rho[i][j]*(1.0 +3.0*tmp +9.0/2.0*tmp*tmp -3.0/2.0*u2)/36.0;
      }
    } }

    for(i = 0; i <= nx; i++){ for(j = 0; j <= ny; j++){ for(k = 0; k <= 8; k++){
      f[k][i][j] = f[k][i][j] - (f[k][i][j] - f0[k][i][j])/tau;
    } } }

    // 一様体積力を x 方向流れの駆動項として加える。
    for(i = 0; i <= nx; i++){ for(j = 0; j <= ny; j++){
      for(k = 1; k <= 4; k++){
        f[k][i][j] = f[k][i][j] + rho[i][j]*(cx[k]*gx + cy[k]*gy)/3.0;
      }
      for(k = 5; k <= 8; k++){
        f[k][i][j] = f[k][i][j] + rho[i][j]*(cx[k]*gx + cy[k]*gy)/12.0;
      }
    } }

    // 衝突後の分布関数を離散速度方向へ streaming する。
    for(i = 0; i <= nx; i++){ for(j = 0; j <= ny; j++){ for(k = 0; k <= 8; k++){
      ftmp[k][i][j] = f[k][i][j];
    } } }

    for(i = 0; i <= nx; i++){ for(j = 0; j <= ny; j++){ for(k = 0; k <= 8; k++){
      in = i + (int)cx[k]; jn = j + (int)cy[k];
      if(in == nx + 1){in =  0;} if(in == -1){in = nx;}
      if(jn == ny + 1){jn =  0;} if(jn == -1){jn = ny;}
      f[k][in][jn] = ftmp[k][i][j];
    } } }

    // 上下壁では bounce-back により no-slip 条件を近似する。
    for(i = 0; i <= nx; i++){
      f[2][i][0] = f[4][i][0]; f[4][i][ny] = f[2][i][ny];
      f[5][i][0] = f[7][i][0]; f[7][i][ny] = f[5][i][ny];
      f[6][i][0] = f[8][i][0]; f[8][i][ny] = f[6][i][ny];
    } 

    // 分布関数から密度と速度を再構成する。
    for(i = 0; i <= nx; i++){ for(j = 0; j <= ny; j++){
      rho[i][j] = f[0][i][j]; u[i][j] = 0; v[i][j] = 0;
      for( k = 1; k <= 8; k++){
        rho[i][j] = rho[i][j] + f[k][i][j];
          u[i][j] =   u[i][j] + f[k][i][j]*cx[k];
          v[i][j] =   v[i][j] + f[k][i][j]*cy[k];
      } 
      u[i][j] = u[i][j]/rho[i][j];
      v[i][j] = v[i][j]/rho[i][j];
    } }

    // 前ステップとの差の最大値を、定常到達の指標として使う。
    norm = 0.0;
    for(i = 0; i <= nx; i++){ for(j = 0; j <= ny; j++){
      tmp = sqrt(pow(u[i][j]-un[i][j],2) + pow(v[i][j]-vn[i][j],2));
      if(tmp > norm){ norm = tmp; }
    } }

  } //loop2

  // 収束の様子、最大流速、壁面 slip を標準出力へ表示する。
  printf("Time = %d, Norm = %15.8e\n", time, norm);
  printf("Umax = %8.6e(%6.4e)\n",u[nx/2][ny/2], gx/8.0/nu*ny*ny);
  printf("Slip = %8.6e(%6.4e)\n",u[nx/2][0]/(gx/8.0/nu*ny*ny),
                                 8.0*tau*(2.0*tau -1.0)/3.0/ny/ny); 

  // 中心断面速度を簡易な横棒グラフとして文字出力する。
  for(i = 0; i <= nx; i++){ for(j = 0; j <= ny; j++){
    a[i][j]='0';
  } }

  for(j = 0; j <= ny; j++){
    nb = u[nx/2][j]/(gx/8.0/nu*h*h)*nx;

    for(i = 0; i <= nb; i++){
      a[i][j]= '-';
    }
  }

  for(j = 0; j <= ny; j++){ for(i = 0; i <= nx; i++){
    printf("%c",a[i][j]);
    }
    printf("\n");
  }

  // 十分収束したら、中心断面の無次元速度分布を data に保存して終了する。
  if(norm < 0.0000000001 && time > 10000){
    fp = fopen("data","w");

    for(j = 0; j <= ny; j++){
      fprintf(fp,"%10.8e\n", u[nx/2][j]/(gx/8.0/nu*h*h));
    }
    fclose(fp);

    exit(0);
  }

  } //loop1

  return 0;
}
