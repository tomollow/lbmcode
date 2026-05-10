// Advection equation 
//
// de/dt + c*de/dx = 0
//
// Finite Difference Method.
// The source code is written in C programming language.
// Copyright @2021 Takeshi Seta All Rights Reserved.
// e  : numerical solution (time = n + 1) 
// en : numerical solution (time = n    ) 
// eo : numerical solution (time = n - 1) 
// dx : grid spacing
// dt : time step
// nx : number of grid points 
// c  : advection velocity
// flag == 1 (Upwind Scheme)
// flag == 2 (FTCS Scheme)
// flag == 3 (Lax-Wendroff Scheme)
// flag == 4 (Leap-Frog Scheme)

#include<stdio.h>
#include<stdlib.h>
#define DIM 100
int main(int argc, char *argv[])
{
  FILE *fp;
  FILE *history_fp;
  int i, j, k, m, n, nx, flag;
  double time, dt, dx, c, e[DIM], en[DIM], eo[DIM];
  char a[50][50];

  dt = 0.5;
// dt = 1.1;
  nx = 50;
  time = 0.0; dx = 1.0; c = 1.0;

  flag = 1;
  if(argc >= 2){
    flag = atoi(argv[1]);
  }else{
    printf("flag?");
    scanf("%d",&flag);
  }
  if(argc >= 3){
    dt = atof(argv[2]);
  }
  if(flag < 1 || flag > 4){
    fprintf(stderr,"flag must be 1 (Upwind), 2 (FTCS), 3 (Lax-Wendroff), or 4 (Leap-Frog).\n");
    return 1;
  }

// initial condition
  for(i = 0; i < nx; i++){
    e[i] = 0.0;
   en[i] = e[i];
   eo[i] = e[i];
  }
  for(i = 1; i <= 5; i++){
    e[i] = 1.0;
   en[i] = e[i];
   eo[i] = e[i];
  }

//calculation
  history_fp = fopen("fdmadv_history.dat","w");
  for(k = 0; k < 6; k++){
    for(m = 0; m < 10; m++){
      time =  time + dt;

      if(flag == 1){ // Upwind
        for(i = 1; i < nx; i++){
          e[i] = en[i] - c*dt/dx*(en[i] - en[i-1]);
        }
      }else if(flag == 2){ // FTCS
        for(i = 1; i < nx - 1; i++){
          e[i] = en[i] - 0.5*c*dt/dx*(en[i+1] - en[i-1]);
        }
      }else if(flag == 3){ // Lax Wendroff
        for(i = 1; i < nx - 1; i++){
          e[i] = en[i] - 0.5*c*dt/dx*(en[i+1] - en[i-1])
          + 0.5*c*c*dt*dt/dx/dx*(en[i+1] - 2.0*en[i] + en[i-1]);
        }
      }else if(flag == 4){ // Leap-Frog
        for(i = 1; i < nx - 1; i++){
          e[i] = eo[i] - c*dt/dx*(en[i+1] - en[i-1]);
        }
      }

      for(i = 0; i < nx; i++){
        eo[i] = en[i];
        en[i] =  e[i];
      }
    }

    for(i = 0; i < nx; i++){
      n = (int)(e[i]*5.0);
      if(n < 0){
        n = 0;
      }
      if(n > 10){
        n = 10;
      }
      for(j = 0; j < 10; j++){
        a[i][j] = ' ';
      }
      for(j = 0; j < n; j++){
        a[i][j] = 'o';
      }
    }
    printf("time = %f\n",time);
    for(j = 10- 1; j >= 0; j--){
      for(i = 0; i < nx; i++){
        printf("%c",a[i][j]);
      }
      printf("\n");
    }
    printf("------------------\n");

    fprintf(history_fp,"%10.8e", time);
    for(i = 0; i < nx; i++){
      fprintf(history_fp," %10.8e", e[i]);
    }
    fprintf(history_fp,"\n");
  }

  fclose(history_fp);

  fp = fopen("fdmadv","w");
  for(i = 0; i < nx; i++){
    fprintf(fp," %10.8e", e[i]);
  } 
  fclose(fp);

  return 0;

}
