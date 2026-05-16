#include <stdio.h>
#include <math.h>
#include "nr.h"
#include "nrutil.h"
#include "nrutil.c"
#include "jacobi.c"
#include "lubksb.c"
#include "ludcmp.c"
#include <time.h>


double *initm0(int n) {
    int i;
    double *m0;
    m0 = dvector(1, n);
    for(i = 1; i <= n; i++)
        m0[i] = /*0.1 + 0.1 * (i - 1)**/ rand() / (double)RAND_MAX;
}

double *initk0(int n) {
    int i;
    double *k0;
    k0 = dvector(1, n);
    for(i = 1; i <= n; i++)
        k0[i] =/* 10.0* i**/100*rand() / (double)RAND_MAX;
}


double *initx0(int n, double min_val, double max_val) {
    srand(time(NULL));
    double *x0 = dvector(1, n);
    for(int i=1; i<=n; i++)
        x0[i] = min_val + (max_val - min_val) * rand() / (double)RAND_MAX;
    return x0;
}

double *initxp0(int n, double min_val, double max_val) {
    srand(time(NULL) + 1); 
    double *xp0 = dvector(1, n);
    for(int i=1; i<=n; i++)
        xp0[i] = min_val + (max_val - min_val) * rand() / (double)RAND_MAX;
    return xp0;
}



double **initk(double *k0, int n) {
    int i, j;
    double **k;
    k = dmatrix(1, n, 1, n);
    
    for (i = 1; i <= n; i++)
        for (j = 1; j <= n; j++)
            if (i == j)
                k[i][j] = k0[i] + (i == n ? 0 : k0[i + 1]);
            else if (i - j == 1)
                k[i][j] = -k0[i];
            else if (j - i == 1)
                k[i][j] = -k0[j];
            else
                k[i][j] = 0;
    return k;
}


double **calck1(double **k, double *m0, int n) {
    int i, j;
    double **k1;
    k1 = dmatrix(1, n, 1, n);
    
    for (i = 1; i <= n; i++)
        for (j = 1; j <= n; j++)
            if (i == j)
                k1[i][j] = k[i][j] / m0[i];
            else if (i - j == 1 || j - i == 1)
                k1[i][j] = k[i][j] / sqrt(m0[i] * m0[j]);
            else
                k1[i][j] = 0;
    return k1;
}


void calcw(double *w, int n) {
    int i;
    for (i = 1; i <= n; i++)
        w[i] = sqrt(w[i]);
}


void calca(double **a, double *m0, int n) {
    int i, j;
    for (i = 1; i <= n; i++)
        for (j = 1; j <= n; j++)
            a[i][j] /= sqrt(m0[i]);
}


void copy(double **to, double **from, int n) {
    int i, j;
    for (i = 1; i <= n; i++)
        for (j = 1; j <= n; j++)
            to[i][j] = from[i][j];
}

void calcbeta(double *xp0, double *w, int n) {
    int i;
    for (i = 1; i <= n; i++)
        xp0[i] /= w[i];
}

int main(void) {
    int i, j, n, nrot, *indx;
    double **a, **alu, *m0, **k, *k0, **k1, *x0, *xp0, *w, d, t;
    FILE *out;
    
    n = 10;  // Количество масс в системе
    
  
    m0 = initm0(n);
    k0 = initk0(n);
    x0 = initx0(n, -0.3, 0.3);
	xp0 = initxp0(n, -0.5, 0.5);
    
   
    k = initk(k0, n);
    k1 = calck1(k, m0, n);

    w = dvector(1, n);
    a = dmatrix(1, n, 1, n);
    
    // с/з
    jacobi(k1, n, w, a, &nrot);
    

    calcw(w, n);
    calca(a, m0, n);
    
    // LU-декомпозиция для решения систем уравнений
    indx = ivector(1, n);
    alu = dmatrix(1, n, 1, n);
    copy(alu, a, n);
    
    ludcmp(alu, n, indx, &d);
    lubksb(alu, n, indx, x0);
    lubksb(alu, n, indx, xp0);
    
    calcbeta(xp0, w, n);
    
    
    out = fopen("osc.txt", "w");
    
    // столбцы
    fprintf(out, "%-12s", "Time"); 
    for (i = 1; i <= n; i++)
        fprintf(out, "%-12s%d", "Mass", i);  
    fprintf(out, "\n");
    
    // расчёты
    for (t = 0; t <= 100; t += 0.5) {
        fprintf(out, "%-12.5e", t);  // Вывод времени
        
        
        for (i = 1; i <= n; i++) {
            double x = 0;
          
            for (j = 1; j <= n; j++) {
                x += a[i][j] * (x0[j] * cos(w[j] * t) + xp0[j] * sin(w[j] * t));
            }
            x += 4 * i; 
            fprintf(out, "%-12.5e", x);  
        }
        fprintf(out, "\n");  
    }
    
    fclose(out);
    
    
    free_dvector(m0, 1, n);
    free_dvector(k0, 1, n);
    free_dvector(x0, 1, n);
    free_dvector(xp0, 1, n);
    free_dvector(w, 1, n);
    free_dmatrix(k, 1, n, 1, n);
    free_dmatrix(k1, 1, n, 1, n);
    free_dmatrix(a, 1, n, 1, n);
    free_dmatrix(alu, 1, n, 1, n);
    free_ivector(indx, 1, n);
    
    return 0;
}