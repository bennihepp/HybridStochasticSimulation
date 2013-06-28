#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>

#include "ch_ethz_khammash_ode_lsodar_LsodarDirectSolver.h"
#include "jni_utils.h"

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

struct lsodar_direct_data {
    CALLBACK_F_PTR f_ptr;
    int neq[3];
    double* y;
    double t;
    double tout;
    int itol;
    double rtol[1];
    double atol[1];
    int itask;
    int istate;
    int iopt;
    double* rwork;
    int lrw;
    double* iwork;
    int liw;
    CALLBACK_JAC_PTR jac_ptr;
    int jt;
    CALLBACK_G_PTR g_ptr;
    int ng;
    int* jroot;
    JNIEnv* env;
    jobject ode;
    jobject ef;
    jmethodID vectorFieldMethodID;
    jmethodID eventValuesMethodID;
    jmethodID reportEventMethodID;
    jdouble* x;
    jdouble* xTmp;
    jdouble* xDot;
    jdouble* g;
    jint* eventIndex;
};

void f_direct_bridge(int neq[], double* t, double y[], double ydot[]) {
    void** ptr = (void**)(&neq[1]);
    struct lsodar_direct_data* data = (struct lsodar_direct_data*)ptr[0];
    JNIEnv *env = data->env;

    // Copy y to the JVM direct array xTmp, then compute the vectorfield and copy
    // the JVM direct array xDot to ydot
    memcpy(data->xTmp, y, neq[0] * sizeof(double));
    (*env)->CallVoidMethod(env, data->ode, data->vectorFieldMethodID, *t);
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "f_direct_bridge: Failed to call computeVectorField!\n");
        return;
    }
    memcpy(ydot, data->xDot, neq[0] * sizeof(double));
}

void g_direct_bridge(int neq[], double* t, double y[], int* ng, double gout[]) {
    void** ptr = (void**)(&neq[1]);
    struct lsodar_direct_data* data = (struct lsodar_direct_data*)ptr[0];
    JNIEnv *env = data->env;

    // Copy y to the JVM direct array xTmp, then compute the event functions
    // and copy the JVM direct array g to gout
    memcpy(data->xTmp, y, neq[0] * sizeof(double));
    (*env)->CallVoidMethod(env, data->ef, data->eventValuesMethodID, *t);
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "g_direct_bridge: Failed to call computeEventValues!\n");
        return;
    }
    memcpy(gout, data->g, ng[0] * sizeof(double));
}

JNIEXPORT jlong JNICALL Java_ch_ethz_khammash_ode_lsodar_LsodarDirectSolver_jni_1initialize
  (JNIEnv *env, jobject obj, jobject ode, jobject ef,
  jobject xBuffer, jobject xTmpBuffer, jobject xDotBuffer, jobject gBuffer, jobject eventIndexBuffer,
  jdouble relTol, jdouble absTol) {
    // Get references to classes and some methods and find out the
    // values of neq (number of equations) and ng (number of event functions)
    jclass odeCls = (*env)->GetObjectClass(env, ode);
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to get class of ode object\n");
        return 0;
    }
    jclass efCls = (*env)->GetObjectClass(env, ef);
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to get class of ef object\n");
        return 0;
    }
    jmethodID methodID = (*env)->GetMethodID(env, odeCls, "getDimensionOfVectorField", "()I");
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to get method ID of getDimensionOfVectorField\n");
        return 0;
    }
    jint neq = (*env)->CallIntMethod(env, ode, methodID);
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to call getDimensionOfVectorField\n");
        return 0;
    }
    methodID = (*env)->GetMethodID(env, efCls, "getNumberOfEventValues", "()I");
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to get method ID of getNumberOfEventValues\n");
        return 0;
    }
    jint ng = (*env)->CallIntMethod(env, ef, methodID);
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to call getNumberOfEventValues\n");
        return 0;
    }
    // Allocate and initialize lsodar_direct_data (see documentation of lsodar)
    struct lsodar_direct_data* data = jni_malloc(env, sizeof(struct lsodar_direct_data));
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to allocate lsodar_direct_data with jni_malloc\n");
        return 0;
    }
    data->f_ptr = &f_direct_bridge;
    data->neq[0] = neq;
    void** ptr = (void**)(&data->neq[1]);
    ptr[0] = data;
    data->y = jni_malloc(env, neq * sizeof(double));
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to allocate y with jni_malloc\n");
        jni_free(env, data);
        return 0;
    }
    data->t = 0;
    data->tout = 1;
    data->itol = 1;
    data->rtol[1] = relTol;
    data->atol[1] = absTol;
    data->itask = 1;
    data->istate = 1;
    data->iopt = 0;
    int lrn = 20 + 16*neq + 3*ng;
    int lrs = 22 + 9*neq + neq*neq + 3*ng;
    data->lrw = max(lrn, lrs);
    data->rwork = jni_malloc(env, data->lrw * sizeof(double));
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to allocate rwork with jni_malloc\n");
        jni_free(env, data);
        return 0;
    }
    data->liw = 20 + neq;
    data->iwork = jni_malloc(env, data->liw * sizeof(double));
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to allocate iwork with jni_malloc\n");
        jni_free(env, data);
        return 0;
    }
    data->jac_ptr = NULL;
    data->jt = 2;
    data->g_ptr = &g_direct_bridge;
    data->ng = ng;
    data->jroot = jni_malloc(env, ng * sizeof(double));
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to allocate jroot with jni_malloc\n");
        jni_free(env, data);
        return 0;
    }
    data->env = env;
    // Acquire method IDs of callbacks and put them into lsodar_direct_data
    data->ode = (*env)->NewGlobalRef(env, ode);
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to create global reference to ode object\n");
        jni_free(env, data);
        return 0;
    }
    data->ef = (*env)->NewGlobalRef(env, ef);
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to create global reference to ef object\n");
        jni_free(env, data);
        return 0;
    }
    data->vectorFieldMethodID = (*env)->GetMethodID(env, odeCls, "computeVectorField", "(D)V");
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to get method ID of computeVectorField\n");
        jni_free(env, data);
        return 0;
    }
    data->eventValuesMethodID = (*env)->GetMethodID(env, efCls, "computeEventValues", "(D)V");
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to get method ID of computeEventValues\n");
        jni_free(env, data);
        return 0;
    }
    // Check capacities of direct buffer objects
    if ((*env)->GetDirectBufferCapacity(env, xBuffer) < neq) {
        fprintf(stderr, "jni_initialize: Direct buffer xBuffer is not big enough\n");
        jni_free(env, data);
        return 0;
    }
    if ((*env)->GetDirectBufferCapacity(env, xTmpBuffer) < neq) {
        fprintf(stderr, "jni_initialize: Direct buffer xTmpBuffer is not big enough\n");
        jni_free(env, data);
        return 0;
    }
    if ((*env)->GetDirectBufferCapacity(env, xDotBuffer) < neq) {
        fprintf(stderr, "jni_initialize: Direct buffer xDotBuffer is not big enough\n");
        jni_free(env, data);
        return 0;
    }
    if ((*env)->GetDirectBufferCapacity(env, gBuffer) < ng) {
        fprintf(stderr, "jni_initialize: Direct buffer gBuffer is not big enough\n");
        jni_free(env, data);
        return 0;
    }
    if ((*env)->GetDirectBufferCapacity(env, eventIndexBuffer) < 1) {
        fprintf(stderr, "jni_initialize: Direct buffer eventIndexBuffer is not big enough\n");
        jni_free(env, data);
        return 0;
    }
    // Get addresses of direct buffer objects
    data->x = (*env)->GetDirectBufferAddress(env, xBuffer);
    //int i;
    //for (i=0; i < neq; i++)
        //data->x[i] = 0.0;
    if (data->x == NULL) {
        fprintf(stderr, "jni_initialize: Failed to get direct address of xBuffer\n");
        jni_free(env, data);
        return 0;
    }
    data->xTmp = (*env)->GetDirectBufferAddress(env, xTmpBuffer);
    if (data->xTmp == NULL) {
        fprintf(stderr, "jni_initialize: Failed to get direct address of xTmpBuffer\n");
        jni_free(env, data);
        return 0;
    }
    data->xDot = (*env)->GetDirectBufferAddress(env, xDotBuffer);
    if (data->xDot == NULL) {
        fprintf(stderr, "jni_initialize: Failed to get direct address of xDotBuffer\n");
        jni_free(env, data);
        return 0;
    }
    data->g = (*env)->GetDirectBufferAddress(env, gBuffer);
    if (data->g == NULL) {
        fprintf(stderr, "jni_initialize: Failed to get direct address of gBuffer\n");
        jni_free(env, data);
        return 0;
    }
    data->eventIndex = (*env)->GetDirectBufferAddress(env, eventIndexBuffer);
    if (data->eventIndex == NULL) {
        fprintf(stderr, "jni_initialize: Failed to get direct address of eventIndexBuffer\n");
        jni_free(env, data);
        return 0;
    }
    // Pass pointer-address of lsodar_direct_data back to Java
    jlong jni_pointer = (jlong)data;
    return jni_pointer;
}

JNIEXPORT void JNICALL Java_ch_ethz_khammash_ode_lsodar_LsodarDirectSolver_jni_1dispose
  (JNIEnv *env, jobject obj, jlong jni_pointer) {
    struct lsodar_direct_data* data = (struct lsodar_direct_data*)jni_pointer;
    jni_free(env, data->y);
    jni_free(env, data->rwork);
    jni_free(env, data->iwork);
    jni_free(env, data->jroot);
    (*env)->DeleteGlobalRef(env, data->ode);
    (*env)->DeleteGlobalRef(env, data->ef);
    jni_free(env, data);
}

JNIEXPORT jdouble JNICALL Java_ch_ethz_khammash_ode_lsodar_LsodarDirectSolver_jni_1integrate
  (JNIEnv *env, jobject obj, jlong jni_pointer, jdouble t0, jdouble t1) {
    struct lsodar_direct_data* data = (struct lsodar_direct_data*)jni_pointer;
    // Reset event indicator
    //data->eventIndex = -1;
    // Copy start and end time for integration
    data->t = t0;
    data->tout = t1;
    // Call lsodar solver
    dlsodar_(
        data->f_ptr, data->neq, data->x, &data->t, &data->tout,
        &data->itol, data->rtol, data->atol, &data->itask,
        &data->istate, &data->iopt, data->rwork, &data->lrw,
        data->iwork, &data->liw, data->jac_ptr, &data->jt,
        data->g_ptr, &data->ng, data->jroot
    );
    if (data->istate < 0) {
        if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
            fprintf(stderr, "jni_integrate: Java exception was thrown during ODE integration\n");
            return NAN;
        }
        // dlsodar failed so we throw an exception
        jclass excpClass = (*env)->FindClass(env, "ch/ethz/khammash/ode/lsodar/Exception");
        if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
            fprintf(stderr, "jni_integrate: Failed to find class ch.ethz.khammash.ode.lsodar.Exception\n");
            return NAN;
        }
        jmethodID constructorMethodID = (*env)->GetMethodID(env, excpClass, "<init>", "(ID)V");
        if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
            fprintf(stderr, "jni_integrate: Failed to get method ID of constructor for ch.ethz.khammash.ode.lsodar.Exception\n");
            // Reset istate for next integration
            data->istate = 1;
            return 0;
        }
        jobject excp = (*env)->NewObject(env, excpClass, constructorMethodID, data->istate, data->t);
        if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
            fprintf(stderr, "jni_integrate: Failed to create ch.ethz.khammash.ode.lsodar.Exception object\n");
            // Reset istate for next integration
            data->istate = 1;
            return NAN;
        }
        int result = (*env)->Throw(env, excp);
        if (result != 0) {
            fprintf(stderr, "jni_integrate: Failed to throw exception\n");
            // Reset istate for next integration
            data->istate = 1;
            return NAN;
        }
        // Reset istate for next integration
        data->istate = 1;
        return NAN;
    }

    if (data->istate == 3) {
        // Report one event to the solver
        int i;
        for (i=0; i < data->ng; i++)
            if (data->jroot[i] == 1) {
                data->eventIndex[0] = i;
                break;
            }
        // Reset istate for next integration
        data->istate = 1;
    }
    return data->t;
}
