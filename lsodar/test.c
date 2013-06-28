#include <stdlib.h>
#include <stdio.h>

#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
	#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

typedef void (*CALLBACK_F_PTR)(int neq[], double* t, double y[], double ydot[]);
typedef void (*CALLBACK_G_PTR)(int neq[], double* t, double y[], int* ng, double gout[]);
typedef void (*CALLBACK_JAC_PTR)(int neq[], double* t, double y[], int* ml, int* mu, double pd[], int* nrowpd);

extern void dlsodar_(
    CALLBACK_F_PTR f, int neq[], double y[], double* t, double* tout, int* itol,
    double rtol[], double atol[], int* itask, int* istate, int* iopt,
    double rwork[], int* lrw, double iwork[], int* liw,
    CALLBACK_JAC_PTR jac, int* jt, CALLBACK_G_PTR g, int* ng, int jroot[]);

void f(int neq[], double* t, double y[], double ydot[]) {
    ydot[0] = 1.0;
}

void g(int neq[], double* t, double y[], int* ng, double gout[]) {
    gout[0] = 0.5 - y[0];
}

int main(int argc, char** argv) {
    CALLBACK_F_PTR f_ptr = &f;
    int neq = 1;
    int ng = 1;
    double y[1] = {0};
    double t = 0;
    double tout = 1;
    int itol = 1;
    double rtol[1] = {1e-6};
    double atol[1] = {1e-3};
    int itask = 1;
    int istate = 1;
    int iopt = 0;
    int lrn = 20 + 16*neq + 3*ng;
    int lrs = 22 + 9*neq + neq*neq + 3*ng;
    int lrw = max(lrn, lrs);
    double* rwork = calloc(lrw, sizeof(double));
    int liw = 20 + neq;
    double* iwork = calloc(liw, sizeof(double));
    CALLBACK_JAC_PTR jac_ptr = NULL;
    int jt = 2;
    CALLBACK_G_PTR g_ptr = &g;
    int* jroot = calloc(ng, sizeof(double));
    printf("y(0s) = %f\n", y[0]);
    dlsodar_(&f, &neq, y, &t, &tout, &itol, rtol, atol, &itask, &istate, &iopt, rwork, &lrw, iwork, &liw, jac_ptr, &jt, g_ptr, &ng, jroot);
    printf("y(1s) = %f\n", y[0]);
    printf("jroot[0] = %d\n", jroot[0]);
    return 0;
}
