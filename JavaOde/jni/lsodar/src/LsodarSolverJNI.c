#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>

#include "ch_ethz_bhepp_ode_lsodar_LsodarSolver.h"
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

void _gfortran_st_write() {
    printf("gfortran message");
    fflush(stdout);
}
void _gfortran_transfer_character_write() {
}
void _gfortran_st_write_done() {
}

extern void dlsodar_(
    CALLBACK_F_PTR f, int neq[], double y[], double* t, double* tout, int* itol,
    double rtol[], double atol[], int* itask, int* istate, int* iopt,
    double rwork[], int* lrw, double iwork[], int* liw,
    CALLBACK_JAC_PTR jac, int* jt, CALLBACK_G_PTR g, int* ng, int jroot[]);

struct lsodar_data {
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
    jdoubleArray jy;
    jdoubleArray jydot;
    jdoubleArray jg;
};

void f_bridge(int neq[], double* t, double y[], double ydot[]) {
    void** ptr = (void**)(&neq[1]);
    struct lsodar_data* data = (struct lsodar_data*)ptr[0];
    JNIEnv *env = data->env;

    jdoubleArray jy = data->jy;
    (*env)->SetDoubleArrayRegion(env, jy, 0, neq[0], y);
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "f_bridge: Failed to set array region!\n");
        return;
    }
    jdoubleArray jydot = data->jydot;
    (*env)->CallVoidMethod(env, data->ode, data->vectorFieldMethodID, *t, jy, jydot);
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "f_bridge: Failed to call computeVectorField!\n");
        return;
    }
    //jdouble *ydottmp = (*env)->GetDoubleArrayElements(env, jydot, 0);
    jdouble *ydottmp = (*env)->GetPrimitiveArrayCritical(env, jydot, 0);
    if (ydottmp == NULL) {
        fprintf(stderr, "f_bridge: Failed to get primitive array of jydot!\n");
        return;
    }
    memcpy(ydot, ydottmp, neq[0] * sizeof(double));
    (*env)->ReleasePrimitiveArrayCritical(env, jydot, ydottmp, JNI_ABORT);
    //(*env)->ReleaseDoubleArrayElements(env, jydot, ydottmp, JNI_ABORT);
}

void g_bridge(int neq[], double* t, double y[], int* ng, double gout[]) {
    void** ptr = (void**)(&neq[1]);
    struct lsodar_data* data = (struct lsodar_data*)ptr[0];
    JNIEnv *env = data->env;

    jdoubleArray jy = data->jy;
    (*env)->SetDoubleArrayRegion(env, jy, 0, data->neq[0], y);
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "g_bridge: Failed to set array region\n");
        return;
    }
    jdoubleArray jg = data->jg;
    (*env)->CallVoidMethod(env, data->ef, data->eventValuesMethodID, *t, jy, jg);
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "g_bridge: Failed to call computeEventValues!\n");
        return;
    }
    //jdouble *gtmp = (*env)->GetDoubleArrayElements(env, jg, 0);
    jdouble *gtmp = (*env)->GetPrimitiveArrayCritical(env, jg, 0);
    if (gtmp == NULL) {
        fprintf(stderr, "g_bridge: Failed to get primitive array of jg!\n");
        return;
    }
    memcpy(gout, gtmp, ng[0] * sizeof(double));
    (*env)->ReleasePrimitiveArrayCritical(env, jg, gtmp, JNI_ABORT);
    //(*env)->ReleaseDoubleArrayElements(env, jg, gtmp, JNI_ABORT);
}

/*
 * Class:     LsodarOdeSolver
 * Method:    jni_initialize
 * Signature: (LLsodarOde;LLsodarTimepointIterator;)V
 */
JNIEXPORT jlong JNICALL Java_ch_ethz_bhepp_ode_lsodar_LsodarSolver_jni_1initialize
  (JNIEnv *env, jobject obj, jobject ode, jobject ef, jdouble relTol, jdouble absTol) {
    // Get references to classes and some methods and find out the
    // values of neq (number of equations) and ng (number of event functions)
    jclass objCls = (*env)->GetObjectClass(env, obj);
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to get class of obj object\n");
        return 0;
    }
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
    // Allocate and initialize lsodar_data (see documentation of lsodar)
    struct lsodar_data* data = jni_malloc(env, sizeof(struct lsodar_data));
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to allocate lsodar_data with jni_malloc\n");
        return 0;
    }
    data->f_ptr = &f_bridge;
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
    data->g_ptr = &g_bridge;
    data->ng = ng;
    data->jroot = jni_malloc(env, ng * sizeof(double));
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to allocate jroot with jni_malloc\n");
        jni_free(env, data);
        return 0;
    }
    data->env = env;
    // Acquire method IDs of callbacks and put them into lsodar_data
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
    data->vectorFieldMethodID = (*env)->GetMethodID(env, odeCls, "computeVectorField", "(D[D[D)V");
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to get method ID of computeVectorField\n");
        jni_free(env, data);
        return 0;
    }
    data->eventValuesMethodID = (*env)->GetMethodID(env, efCls, "computeEventValues", "(D[D[D)V");
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to get method ID of computeEventValues\n");
        jni_free(env, data);
        return 0;
    }
    data->reportEventMethodID = (*env)->GetMethodID(env, objCls, "reportEvent", "(ID[D)V");
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to get method ID of reportEvent\n");
        jni_free(env, data);
        return 0;
    }
    // Create some java arrays for copying into and out of the JVM
    jdoubleArray jy = (*env)->NewDoubleArray(env, neq);
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to create jy array\n");
        return 0;
    }
    data->jy = (*env)->NewGlobalRef(env, jy);
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to create global reference to jy array\n");
        return 0;
    }
    jdoubleArray jydot = (*env)->NewDoubleArray(env, neq);
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to create jydot array\n");
        return 0;
    }
    data->jydot = (*env)->NewGlobalRef(env, jydot);
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to create global reference to jydot array\n");
        return 0;
    }
    jdoubleArray jg = (*env)->NewDoubleArray(env, ng);
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to create jg array\n");
        return 0;
    }
    data->jg = (*env)->NewGlobalRef(env, jg);
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to create global reference to jg array\n");
        return 0;
    }
    // Pass pointer-address of lsodar_data back to Java
    jlong jni_pointer = (jlong)data;
    return jni_pointer;
}

/*
 * Class:     LsodarOdeSolver
 * Method:    jni_dispose
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_ch_ethz_bhepp_ode_lsodar_LsodarSolver_jni_1dispose
  (JNIEnv *env, jobject obj, jlong jni_pointer) {
    struct lsodar_data* data = (struct lsodar_data*)jni_pointer;
    jni_free(env, data->y);
    jni_free(env, data->rwork);
    jni_free(env, data->iwork);
    jni_free(env, data->jroot);
    (*env)->DeleteGlobalRef(env, data->ode);
    (*env)->DeleteGlobalRef(env, data->ef);
    /*(*env)->DeleteGlobalRef(env, data->jy);
    (*env)->DeleteGlobalRef(env, data->jydot);
    (*env)->DeleteGlobalRef(env, data->jg);*/
    jni_free(env, data);
}

/*
 * Class:     LsodarOdeSolver
 * Method:    jni_integrate
 * Signature: (D[DD)D
 */
JNIEXPORT jdouble JNICALL Java_ch_ethz_bhepp_ode_lsodar_LsodarSolver_jni_1integrate
  (JNIEnv *env, jobject obj, jlong jni_pointer, jdouble t0, jobject jy0, jdouble t1) {
    struct lsodar_data* data = (struct lsodar_data*)jni_pointer;
    // Copy initial state vector
    //jdouble *y0 = (*env)->GetDoubleArrayElements(env, jy0, 0);
    jdouble *y0 = (*env)->GetPrimitiveArrayCritical(env, jy0, 0);
    if (y0 == NULL) {
        fprintf(stderr, "jni_initialize: Failed to get primitive array of jy0\n");
        return NAN;
    }
    memcpy(data->y, y0, data->neq[0] * sizeof(double));
    //(*env)->ReleasePrimitiveArrayCritical(env, jy0, y0, JNI_ABORT);
    //(*env)->ReleaseDoubleArrayElements(env, jy0, y0, JNI_ABORT);

    // Copy start and end time for integration
    data->t = t0;
    data->tout = t1;
    // Call lsodar solver
    dlsodar_(
        data->f_ptr, data->neq, data->y, &data->t, &data->tout,
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
        jclass excpClass = (*env)->FindClass(env, "ch/ethz/bhepp/ode/lsodar$JniException");
        if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
            fprintf(stderr, "jni_integrate: Failed to find class ch.ethz.bhepp.ode.lsodar.JniException\n");
            return NAN;
        }
        jmethodID constructorMethodID = (*env)->GetMethodID(env, excpClass, "<init>", "(ID)V");
        if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
            fprintf(stderr, "jni_integrate: Failed to get method ID of constructor for ch.ethz.bhepp.ode.lsodar.JniException\n");
            // Reset istate for next integration
            data->istate = 1;
            return 0;
        }
        jobject excp = (*env)->NewObject(env, excpClass, constructorMethodID, data->istate, data->t);
        if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
            fprintf(stderr, "jni_integrate: Failed to create ch.ethz.bhepp.ode.lsodar.JniException object\n");
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

    //y0 = (*env)->GetDoubleArrayElements(env, data->jy0, 0);
    //y0 = (*env)->GetPrimitiveArrayCritical(env, jy0, 0);
    //if (y0 == NULL) {
    //    fprintf(stderr, "jni_initialize: Failed to get primitive array of jy0\n");
    //    return NAN;
    //}
    memcpy(y0, data->y, data->neq[0] * sizeof(double));
    (*env)->ReleasePrimitiveArrayCritical(env, jy0, y0, 0);
    //(*env)->ReleaseDoubleArrayElements(env, data->jy, y, 0);

    if (data->istate == 3) {
        // Report event to callback
        int i;
        for (i=0; i < data->ng; i++)
            if (data->jroot[i] == 1) {
                (*env)->CallVoidMethod(env, obj, data->reportEventMethodID, i, data->t, data->jy);
                if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
                    fprintf(stderr, "jni_initialize: Failed to call reportEventMethod\n");
                    // Reset istate for next integration
                    data->istate = 1;
                    return NAN;
                }
            }
        // Reset istate for next integration
        data->istate = 1;
    }
    return data->t;
}
