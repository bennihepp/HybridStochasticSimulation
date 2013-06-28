/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2010/12/01 22:51:32 $
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Example problem:
 * 
 * The following is a simple example problem, with the coding
 * needed for its solution by CVODE. The problem is from
 * chemical kinetics, and consists of the following three rate
 * equations:         
 *    dy1/dt = -.04*y1 + 1.e4*y2*y3
 *    dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*(y2)^2
 *    dy3/dt = 3.e7*(y2)^2
 * on the interval from t = 0.0 to t = 4.e10, with initial
 * conditions: y1 = 1.0, y2 = y3 = 0. The problem is stiff.
 * While integrating the system, we also use the rootfinding
 * feature to find the points at which y1 = 1e-4 or at which
 * y3 = 0.01. This program solves the problem with the BDF method,
 * Newton iteration with the CVDENSE dense linear solver, and a
 * user-supplied Jacobian routine.
 * It uses a user-supplied function to compute the error weights
 * required for the WRMS norm calculations.
 * Output is printed in decades from t = .4 to t = 4.e10.
 * Run statistics (optional outputs) are printed at the end.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <cvode/cvode.h>             /* prototypes for CVODE fcts. and consts. */
#include <cvode/cvode_impl.h>  
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, functions, and macros */
#include <sundials/sundials_dense.h> /* definitions DlsMat and DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */
#include <sundials/sundials_math.h>  /* definition of ABS */

extern int N_SPECIES;
extern int N_RXNS;
extern int N_PARAMS;

/* User-defined vector and matrix accessor macros: Ith, IJth */

/* These macros are defined in order to write code which exactly matches
   the mathematical problem description given above.

   Ith(v,i) references the ith component of the vector v, where i is in
   the range [1..NEQ] and NEQ is defined below. The Ith macro is defined
   using the N_VIth macro in nvector.h. N_VIth numbers the components of
   a vector starting from 0.

   IJth(A,i,j) references the (i,j)th element of the dense matrix A, where
   i and j are in the range [1..NEQ]. The IJth macro is defined using the
   DENSE_ELEM macro in dense.h. DENSE_ELEM numbers rows and columns of a
   dense matrix starting from 0. */

#define Ith(v,i)    NV_Ith_S(v,i)       /* Ith numbers components 1..NEQ */


/* Problem Constants */

#define ZERO   RCONST(0.0)
#define ONE    RCONST(1.0)
#define ETAMX1 RCONST(10000.0) 
#define T0    RCONST(0.0)
#define T1    RCONST(1000.0)
#define NOUT  10000


/* Functions Called by the Solver */

static int fprime(realtype t, N_Vector y, N_Vector ydot, void *user_data);

static int g(realtype t, N_Vector y, realtype *gout, void *user_data);

/* Private functions to output results */

static void PrintOutput(realtype t, realtype y1, realtype y2, realtype y3);
static void PrintRootInfo(int root_f1, int root_f2);

/* Private function to print final statistics */

static void PrintFinalStats(void *cvode_mem);

/* Private function to check function return values */

static int check_flag(void *flagvalue, char *funcname, int opt);


/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

double g_lambda;
double* g_params;
double* g_prop;
int g_deterministic_rxn[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 };
//int g_deterministic_rxn[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15 };
//int g_deterministic_rxn[] = { };
int* g_deltaX;

int main()
{
    int neq = N_SPECIES + 3;
    realtype reltol, abstol;
    realtype t, tout, dt;
    N_Vector y;
    void *cvode_mem;
    int flag, flagr, iout;
    int rootsfound[2];

    t = T0;
    reltol = 1e-6;
    abstol = 1e-6;
    y = NULL;
    cvode_mem = NULL;
    dt = (T1 - T0) / NOUT;

    g_params = malloc(sizeof(double) * N_PARAMS);
    initialize_parameters(g_params);
    g_deltaX = malloc(sizeof(int) * N_SPECIES);
    g_prop = malloc(sizeof(double) * N_RXNS);

    /* Create serial vector of length NEQ for I.C. */
    void* ydata = malloc(sizeof(realtype) * neq);
    y = N_VMake_Serial(neq, ydata);
    if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);

    initialize_state(NV_DATA_S(y));

    double u = rand() / (double)RAND_MAX;
    double lambda = -log(u);
    //Ith(y, N_SPECIES) = 0.0;
    //Ith(y, N_SPECIES+1) = lambda;
    Ith(y, N_SPECIES) = 0.0;
    g_lambda = lambda;

    /* Call CVodeCreate to create the solver memory and specify the 
    * Backward Differentiation Formula and the use of a Newton iteration */
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);

    /* Call CVodeInit to initialize the integrator memory and specify the
    * user's right hand side function in y'=f(t,y), the inital time T0, and
    * the initial dependent variable vector y. */
    flag = CVodeInit(cvode_mem, fprime, t, y);
    if (check_flag(&flag, "CVodeInit", 1)) return(1);

    /* Call CVodeSStolerances to specify the scalar relative tolerance
    * and scalar absolute tolerance */
    flag = CVodeSStolerances(cvode_mem, reltol, abstol);
    if (check_flag(&flag, "CVodeSStolerances", 1)) return(1);

    /* Call CVodeRootInit to specify the root function g with 2 components */
    flag = CVodeRootInit(cvode_mem, 1, g);
    if (check_flag(&flag, "CVodeRootInit", 1)) return(1);

    /* Call CVDense to specify the CVDENSE dense linear solver */
    flag = CVDense(cvode_mem, neq);
    if (check_flag(&flag, "CVDense", 1)) return(1);

    flag = CVodeSetMaxNumSteps(cvode_mem, 1000);
    if (check_flag(&flag, "CVodeSetMaxNumSteps", 1)) return(1);

    /* In loop, call CVode, print results, and test for error.
     Break out of loop when NOUT preset output times have been reached.  */
    printf("\nPDMP\n\n");

    int stochCnt = 0;
    iout = 0;
    tout = t + dt;
    int rootCnt = 0;
    FILE* f = fopen("abc.txt", "w");
    if (f == NULL) {
        perror("Unable to open file");
        return -1;
    }

    srand(time(NULL));

    clock_t start = clock();
    while(t < T1) {

        int i;
        int rootfound = 0;

        //if (sizeof(g_deterministic_rxn) / sizeof(g_deterministic_rxn[0]) > 0) {
            flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
            //PrintOutput(t, Ith(y,1), Ith(y,2), Ith(y,3));

            if (flag != CV_SUCCESS && flag != CV_ROOT_RETURN)
                return 1;

            if (flag == CV_ROOT_RETURN) {
                rootfound = 1;
                flagr = CVodeGetRootInfo(cvode_mem, rootsfound);
                if (flagr != CV_SUCCESS)
                    return 1;
            }
            //if (Ith(y, N_SPECIES+1) - Ith(y, N_SPECIES) < 0)
            //    rootfound = 1;
            tout += dt;
        /*} else {
            rootfound = 1;
            propensity(NV_DATA_S(y), g_params, g_prop);
            for (i=0; i < sizeof(g_deterministic_rxn) / sizeof(g_deterministic_rxn[0]); i++) {
                int rxn = g_deterministic_rxn[i];
                g_prop[rxn] = 0.0;
            }

            double propSum = 0.0;
            for (i=0; i < N_RXNS; i++)
                propSum += g_prop[i];
            u = rand() / (double)RAND_MAX;
            t = tout - log(u) / propSum;
            tout = t;
        }*/

        if (rootfound) {
            
            //PrintRootInfo(rootsfound[0],rootsfound[1]);

            //fprintf(f, "root: %f", t);
            //for (i=0; i < neq; i++)
            //    fprintf(f, "\t%f", Ith(y, i));
            //fprintf(f, "\n");

            int rxn = -1;
            propensity(NV_DATA_S(y), g_params, g_prop);
            for (i=0; i < sizeof(g_deterministic_rxn) / sizeof(g_deterministic_rxn[0]); i++) {
                int rxn = g_deterministic_rxn[i];
                g_prop[rxn] = 0.0;
            }

            double propSum = 0.0;
            for (i=0; i < N_RXNS; i++)
                propSum += g_prop[i];
            u = propSum * rand() / RAND_MAX;
	        double w = 0.0;
	        for (i=0; i < N_RXNS; i++) {
	        	w = w + g_prop[i];
	        	if (u < w) {
	        		rxn = i;
                    break;
	        	}
	        }

            if (rxn == -1)
                return 1;

            CVodeMem cv_mem = (CVodeMem)cvode_mem;

            state_change(g_deltaX, rxn);
            for (i=0; i < N_SPECIES; i++)
              //Ith(y,i) += g_deltaX[i];
              //Ith(cv_mem->cv_zn[0],i) += g_deltaX[i];

            u = rand() / (double)RAND_MAX;
            lambda = -log(u);
            g_lambda = Ith(y, N_SPECIES) + lambda;
            //Ith(y, N_SPECIES) = 0.0;
            //Ith(y, N_SPECIES+1) = lambda;
            //printf("lambda=%f\n", lambda);
            //Ith(y, N_SPECIES+2) += 300.0;

            rootCnt++;
            //flag = CVodeReInit(cvode_mem, t, y) != CV_SUCCESS;
            //if (flag != CV_SUCCESS)
            //    return 1;
            if (cv_mem != NULL) {
                //cv_mem->cv_tn = t;
                N_VScale(ONE, y, cv_mem->cv_zn[0]);
                //N_VScale(cv_mem->cv_h, cv_mem->cv_zn[1], cv_mem->cv_zn[1]);
                /*cv_mem->cv_tn = t;
                cv_mem->cv_q      = 1;
                cv_mem->cv_L      = 2;
                cv_mem->cv_qwait  = cv_mem->cv_L;
                cv_mem->cv_etamax = ETAMX1;
                cv_mem->cv_qu    = 0;
                cv_mem->cv_hu    = ZERO;
                cv_mem->cv_tolsf = ONE;
                N_VScale(ONE, y, cv_mem->cv_zn[0]);
                cv_mem->cv_nst     = 0;
                cv_mem->cv_nfe     = 0;
                cv_mem->cv_ncfn    = 0;
                cv_mem->cv_netf    = 0;
                cv_mem->cv_nni     = 0;
                cv_mem->cv_nsetups = 0;
                cv_mem->cv_nhnil   = 0;
                cv_mem->cv_nstlp   = 0;
                cv_mem->cv_nscon   = 0;
                cv_mem->cv_nge     = 0;
                cv_mem->cv_irfnd   = 0;
                cv_mem->cv_h0u      = ZERO;
                cv_mem->cv_next_h   = ZERO;
                cv_mem->cv_next_q   = 0;
                cv_mem->cv_nor = 0;
                int k;
                for (i = 1; i <= 5; i++)
                    for (k = 1; k <= 3; k++) 
                        cv_mem->cv_ssdat[i-1][k-1] = ZERO;*/
            }

            stochCnt += 1;
  
        } else {
            iout++;
            int i;
            fprintf(f, "out: %f", t);
            for (i=0; i < neq; i++)
                fprintf(f, "\t%f", Ith(y, i));
            fprintf(f, "\n");
        }

    }

    clock_t end = clock();
    fclose(f);
    double runtime = 1000.0 * (end - start) / (double)CLOCKS_PER_SEC;
    printf("Execution time: %f\n", runtime);
    printf("root count: %d\n", rootCnt);
    printf("reaction count: %d\n", stochCnt);

    N_Vector ydot;
    ydot = NULL;
    void* ydotdata = malloc(sizeof(realtype) * neq);
    ydot = N_VMake_Serial(neq, ydotdata);
    if (check_flag((void *)ydot, "N_VNew_Serial", 0)) return(1);

    int i;

    printf("y:");
    for (i=0; i < N_SPECIES; i++)
        printf("\t%f", Ith(y, i));
    printf("\n");

    fprime(t, y, ydot, NULL);
    printf("ydot:");
    for (i=0; i < N_SPECIES; i++)
        printf("\t%f", Ith(ydot, i));
    printf("\n");

    propensity(NV_DATA_S(y), g_params, g_prop);
    printf("prop:");
    for (i=0; i < N_RXNS; i++)
        printf("\t%f", g_prop[i]);
    printf("\n");

    /* Print some final statistics */
    //PrintFinalStats(cvode_mem);

    /* Free y vector */
    N_VDestroy_Serial(y);

    /* Free integrator memory */
    CVodeFree(&cvode_mem);

    return(0);
}


/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */

/*
 * f routine. Compute function f(t,y).
 */

static int fprime(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    int i;
    int j;
    propensity(NV_DATA_S(y), g_params, g_prop);
    for (j=0; j < N_SPECIES; j++)
        Ith(ydot,j) = 0.0;
    for (i=0; i < sizeof(g_deterministic_rxn) / sizeof(g_deterministic_rxn[0]); i++) {
        int rxn = g_deterministic_rxn[i];
        state_change(g_deltaX, rxn);
        for (j=0; j < N_SPECIES; j++) {
            Ith(ydot,j) += g_prop[rxn] * g_deltaX[j];
        }
        g_prop[rxn] = 0.0;
    }
    double propSum = 0.0;
    for (i=0; i < N_RXNS; i++)
        propSum += g_prop[i];

    Ith(ydot,N_SPECIES) = propSum;

    return 0;
}

/*
 * g routine. Compute functions g_i(t,y) for i = 0,1. 
 */

static int g(realtype t, N_Vector y, realtype *gout, void *user_data)
{
    //gout[0] = Ith(y, N_SPECIES+1) - Ith(y, N_SPECIES);
    gout[0] = g_lambda - Ith(y, N_SPECIES);
    //gout[0] = 1.0;
    return 0;
}

/*
 *-------------------------------
 * Private helper functions
 *-------------------------------
 */

static void PrintOutput(realtype t, realtype y1, realtype y2, realtype y3)
{
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("At t = %0.4Le      y =%14.6Le  %14.6Le  %14.6Le\n", t, y1, y2, y3);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("At t = %0.4le      y =%14.6le  %14.6le  %14.6le\n", t, y1, y2, y3);
#else
  printf("At t = %0.4e      y =%14.6e  %14.6e  %14.6e\n", t, y1, y2, y3);
#endif

  return;
}

static void PrintRootInfo(int root_f1, int root_f2)
{
  printf("    rootsfound[] = %3d %3d\n", root_f1, root_f2);

  return;
}

/* 
 * Get and print some final statistics
 */

static void PrintFinalStats(void *cvode_mem)
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  int flag;

  flag = CVodeGetNumSteps(cvode_mem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

  flag = CVDlsGetNumJacEvals(cvode_mem, &nje);
  check_flag(&flag, "CVDlsGetNumJacEvals", 1);
  flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
  check_flag(&flag, "CVDlsGetNumRhsEvals", 1);

  flag = CVodeGetNumGEvals(cvode_mem, &nge);
  check_flag(&flag, "CVodeGetNumGEvals", 1);

  printf("\nFinal Statistics:\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
	 nst, nfe, nsetups, nfeLS, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n",
	 nni, ncfn, netf, nge);
}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer 
 */

static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}
