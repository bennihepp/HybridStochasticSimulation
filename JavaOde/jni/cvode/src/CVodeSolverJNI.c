#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>

#include <cvode/cvode.h>
#include <cvode/cvode_impl.h>
#include <cvode/cvode_dense.h>
#include <nvector/nvector_serial.h>

#include "ch_ethz_bhepp_ode_cvode_CVodeSolver.h"
#include "jni_utils.h"

//typedef int (*CALLBACK_F_PTR)(double t, N_Vector y, N_Vector ydot, void* user_data);
typedef int (*CALLBACK_G_PTR)(double t, N_Vector y, double* gout, void* user_data);

struct cvode_solver_data {
    int neq;
    int ng;
    void* cvode_mem;
    N_Vector y;
    N_Vector yTmp;
    int* rootsfound;
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
    //jint* eventIndex;
    jint eventIndex;
};

int f_bridge(double t, N_Vector y, N_Vector ydot, void* user_data) {
    struct cvode_solver_data* data = (struct cvode_solver_data*)user_data;
    JNIEnv *env = data->env;

    // Copy y to the JVM direct array xTmp, then compute the vectorfield and copy
    // the JVM direct array xDot to ydot
    memcpy(data->xTmp, NV_DATA_S(y), data->neq * sizeof(double));
    (*env)->CallVoidMethod(env, data->ode, data->vectorFieldMethodID, t);
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "f_direct_bridge: Failed to call computeVectorField!\n");
        return -1;
    }
    memcpy(NV_DATA_S(ydot), data->xDot, data->neq * sizeof(double));
    return 0;
}

int g_bridge(double t, N_Vector y, double* gout, void* user_data) {
    struct cvode_solver_data* data = (struct cvode_solver_data*)user_data;

    if (data->ng == 0)
        return 0;

    JNIEnv *env = data->env;

    // Copy y to the JVM direct array xTmp, then compute the event functions
    // and copy the JVM direct array g to gout
    memcpy(data->xTmp, NV_DATA_S(y), data->neq * sizeof(double));
    (*env)->CallVoidMethod(env, data->ef, data->eventValuesMethodID, t);
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "g_direct_bridge: Failed to call computeEventValues!\n");
        return -1;
    }
    memcpy(gout, data->g, data->ng * sizeof(double));
    return 0;
}

/*int jac_bridge(long int N, double t, N Vector y, N Vector fy, DlsMat Jac,
               void *user data, N Vector tmp1, N Vector tmp2, N Vector tmp3) {
    return 0;
}*/

int throw_java_exception(JNIEnv* env, char* msg) {
    jclass excpClass = (*env)->FindClass(env, "ch/ethz/bhepp/ode/cvode/CVodeSolver$JniException");
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "throw_java_exception: Failed to find class ch.ethz.bhepp.ode.cvode.CVodeSolver.JniException\n");
        fprintf(stderr, "Original message: %s\n", msg);
        return -1;
    }
    jmethodID constructorMethodID = (*env)->GetMethodID(env, excpClass, "<init>", "(Ljava/lang/String;)V");
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "throw_java_exception: Failed to get method ID of constructor for ch.ethz.bhepp.ode.cvode.CVodeSolver.JniException\n");
        fprintf(stderr, "Original message: %s\n", msg);
        return -1;
    }
    jobject excp = (*env)->NewObject(env, excpClass, constructorMethodID, msg);
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "throw_java_exception: Failed to create ch.ethz.bhepp.ode.cvode.CVodeSolver.JniException object\n");
        fprintf(stderr, "Original message: %s\n", msg);
        return -1;
    }
    int result = (*env)->Throw(env, excp);
    if (result != 0) {
        fprintf(stderr, "throw_java_exception: Failed to throw exception\n");
        fprintf(stderr, "Original message: %s\n", msg);
        return -1;
    }
    return 0;
}

int throw_java_exception_error_code(JNIEnv* env, char* msg, int error_code) {
    jclass excpClass = (*env)->FindClass(env, "ch/ethz/bhepp/ode/cvode/CVodeSolver$JniException");
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "throw_java_exception_error_code: Failed to find class ch.ethz.bhepp.ode.cvode.CVodeSolver.JniException\n");
        fprintf(stderr, "Original message: %s\n", msg);
        return -1;
    }
    jmethodID constructorMethodID = (*env)->GetMethodID(env, excpClass, "<init>", "(Ljava/lang/String;I)V");
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "throw_java_exception_error_code: Failed to get method ID of constructor for ch.ethz.bhepp.ode.cvode.CVodeSolver.JniException\n");
        fprintf(stderr, "Original message: %s\n", msg);
        return -1;
    }
    jstring jmsg = (*env)->NewStringUTF(env, msg);
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "throw_java_exception_error_code: Failed to convert message to Java String\n");
        fprintf(stderr, "Original message: %s\n", msg);
        return -1;
    }
    jobject excp = (*env)->NewObject(env, excpClass, constructorMethodID, jmsg, error_code);
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "throw_java_exception_error_code: Failed to create ch.ethz.bhepp.ode.cvode.CVodeSolver.JniException object\n");
        fprintf(stderr, "Original message: %s\n", msg);
        return -1;
    }
    int result = (*env)->Throw(env, excp);
    if (result != 0) {
        fprintf(stderr, "throw_java_exception_error_code: Failed to throw exception\n");
        fprintf(stderr, "Original message: %s\n", msg);
        return -1;
    }
    return 0;
}

JNIEXPORT jlong JNICALL Java_ch_ethz_bhepp_ode_cvode_CVodeSolver_jni_1initialize
  (JNIEnv *env, jobject obj, jobject ode, jobject ef,
  jobject xBuffer, jobject xTmpBuffer, jobject xDotBuffer, jobject gBuffer,
  jdouble relTol, jdouble absTol, jint multistepType, jint iterationType,
  jint maxNumOfSteps, jdouble minStep, jdouble maxStep) {
    // Get references to classes and some methods and find out the
    // values of neq (number of equations) and ng (number of event functions)
    jclass odeCls = (*env)->GetObjectClass(env, ode);
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to get class of ode object\n");
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

    jint ng = 0;
    jclass efCls = NULL;
    if (ef != NULL) {
        efCls = (*env)->GetObjectClass(env, ef);
        if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
            fprintf(stderr, "jni_initialize: Failed to get class of ef object\n");
            return 0;
        }
        methodID = (*env)->GetMethodID(env, efCls, "getNumberOfEventValues", "()I");
        if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
            fprintf(stderr, "jni_initialize: Failed to get method ID of getNumberOfEventValues\n");
            return 0;
        }
        ng = (*env)->CallIntMethod(env, ef, methodID);
        if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
            fprintf(stderr, "jni_initialize: Failed to call getNumberOfEventValues\n");
            return 0;
        }
    }

    // Allocate and initialize cvode_solver_data
    struct cvode_solver_data* data = jni_malloc(env, sizeof(struct cvode_solver_data));
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to allocate cvode_solver_data with jni_malloc\n");
        return 0;
    }
    memset(data, 0, sizeof(struct cvode_solver_data));

    data->neq = neq;
    data->ng = ng;

    data->env = env;
    // Acquire method IDs of callbacks and put them into cvode_solver_data
    data->ode = (*env)->NewGlobalRef(env, ode);
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to create global reference to ode object\n");
        jni_free(env, data);
        return 0;
    }
    data->vectorFieldMethodID = (*env)->GetMethodID(env, odeCls, "computeVectorField", "(D)V");
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to get method ID of computeVectorField\n");
        // Delete global java references
        (*env)->DeleteGlobalRef(env, data->ode);
        jni_free(env, data);
        return 0;
    }
    if (ef != NULL) {
        data->ef = (*env)->NewGlobalRef(env, ef);
        if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
            fprintf(stderr, "jni_initialize: Failed to create global reference to ef object\n");
            // Delete global java references
            (*env)->DeleteGlobalRef(env, data->ode);
            jni_free(env, data);
            return 0;
        }
        data->eventValuesMethodID = (*env)->GetMethodID(env, efCls, "computeEventValues", "(D)V");
        if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
            fprintf(stderr, "jni_initialize: Failed to get method ID of computeEventValues\n");
            // Delete global java references
            (*env)->DeleteGlobalRef(env, data->ode);
            (*env)->DeleteGlobalRef(env, data->ef);
            jni_free(env, data);
            return 0;
        }
    } else
        data->ef = NULL;

    // Check capacities of direct buffer objects
    if ((*env)->GetDirectBufferCapacity(env, xBuffer) < neq) {
        throw_java_exception(env, "jni_initialize: Direct buffer xBuffer is not big enough\n");
        // Delete global java references
        (*env)->DeleteGlobalRef(env, data->ode);
        if (data->ef != NULL)
            (*env)->DeleteGlobalRef(env, data->ef);
        jni_free(env, data);
        return 0;
    }
    if ((*env)->GetDirectBufferCapacity(env, xTmpBuffer) < neq) {
        throw_java_exception(env, "jni_initialize: Direct buffer xTmpBuffer is not big enough\n");
        // Delete global java references
        (*env)->DeleteGlobalRef(env, data->ode);
        if (data->ef != NULL)
            (*env)->DeleteGlobalRef(env, data->ef);
        jni_free(env, data);
        return 0;
    }
    if ((*env)->GetDirectBufferCapacity(env, xDotBuffer) < neq) {
        throw_java_exception(env, "jni_initialize: Direct buffer xDotBuffer is not big enough\n");
        // Delete global java references
        (*env)->DeleteGlobalRef(env, data->ode);
        if (data->ef != NULL)
            (*env)->DeleteGlobalRef(env, data->ef);
        jni_free(env, data);
        return 0;
    }
    if ((*env)->GetDirectBufferCapacity(env, gBuffer) < ng) {
        throw_java_exception(env, "jni_initialize: Direct buffer gBuffer is not big enough\n");
        // Delete global java references
        (*env)->DeleteGlobalRef(env, data->ode);
        if (data->ef != NULL)
            (*env)->DeleteGlobalRef(env, data->ef);
        jni_free(env, data);
        return 0;
    }
    /*if ((*env)->GetDirectBufferCapacity(env, eventIndexBuffer) < 1) {
        throw_java_exception(env, "jni_initialize: Direct buffer eventIndexBuffer is not big enough\n");
        // Delete global java references
        (*env)->DeleteGlobalRef(env, data->ode);
        if (data->ef != NULL)
            (*env)->DeleteGlobalRef(env, data->ef);
        jni_free(env, data);
        return 0;
    }*/
    // Get addresses of direct buffer objects
    data->x = (*env)->GetDirectBufferAddress(env, xBuffer);
    if (data->x == NULL) {
        fprintf(stderr, "jni_initialize: Failed to get direct address of xBuffer\n");
        // Delete global java references
        (*env)->DeleteGlobalRef(env, data->ode);
        if (data->ef != NULL)
            (*env)->DeleteGlobalRef(env, data->ef);
        jni_free(env, data);
        return 0;
    }
    data->xTmp = (*env)->GetDirectBufferAddress(env, xTmpBuffer);
    if (data->xTmp == NULL) {
        fprintf(stderr, "jni_initialize: Failed to get direct address of xTmpBuffer\n");
        // Delete global java references
        (*env)->DeleteGlobalRef(env, data->ode);
        if (data->ef != NULL)
            (*env)->DeleteGlobalRef(env, data->ef);
        jni_free(env, data);
        return 0;
    }
    data->xDot = (*env)->GetDirectBufferAddress(env, xDotBuffer);
    if (data->xDot == NULL) {
        fprintf(stderr, "jni_initialize: Failed to get direct address of xDotBuffer\n");
        // Delete global java references
        (*env)->DeleteGlobalRef(env, data->ode);
        if (data->ef != NULL)
            (*env)->DeleteGlobalRef(env, data->ef);
        jni_free(env, data);
        return 0;
    }
    data->g = (*env)->GetDirectBufferAddress(env, gBuffer);
    if (data->g == NULL) {
        fprintf(stderr, "jni_initialize: Failed to get direct address of gBuffer\n");
        // Delete global java references
        (*env)->DeleteGlobalRef(env, data->ode);
        if (data->ef != NULL)
            (*env)->DeleteGlobalRef(env, data->ef);
        jni_free(env, data);
        return 0;
    }
    /*data->eventIndex = (*env)->GetDirectBufferAddress(env, eventIndexBuffer);
    if (data->eventIndex == NULL) {
        fprintf(stderr, "jni_initialize: Failed to get direct address of eventIndexBuffer\n");
        // Delete global java references
        (*env)->DeleteGlobalRef(env, data->ode);
        if (data->ef != NULL)
            (*env)->DeleteGlobalRef(env, data->ef);
        jni_free(env, data);
        return 0;
    }*/

    /* Call CVodeCreate to create the solver memory and specify the 
    * Backward Differentiation Formula and the use of a Newton iteration */
    data->cvode_mem = CVodeCreate(multistepType, iterationType);
    if (data->cvode_mem == NULL) {
        fprintf(stderr, "jni_initialize: Failed to allocate cvode_mem with jni_malloc\n");
        // Delete global java references
        (*env)->DeleteGlobalRef(env, data->ode);
        if (data->ef != NULL)
            (*env)->DeleteGlobalRef(env, data->ef);
        jni_free(env, data);
        return 0;
    }

    if (CVodeSetMaxNumSteps(data->cvode_mem, maxNumOfSteps) != CV_SUCCESS) {
        throw_java_exception(env, "jni_initialize: Failed to set maximum number of steps\n");
        // Free integrator memory
        CVodeFree(&data->cvode_mem);
        // Delete global java references
        (*env)->DeleteGlobalRef(env, data->ode);
        if (data->ef != NULL)
            (*env)->DeleteGlobalRef(env, data->ef);
        jni_free(env, data);
        return 0;
    }

    if (minStep >= 0.0)
        if (CVodeSetMinStep(data->cvode_mem, minStep) != CV_SUCCESS) {
            throw_java_exception(env, "jni_initialize: Failed to set minimum step size\n");
            // Free integrator memory
            CVodeFree(&data->cvode_mem);
            // Delete global java references
            (*env)->DeleteGlobalRef(env, data->ode);
            if (data->ef != NULL)
                (*env)->DeleteGlobalRef(env, data->ef);
            jni_free(env, data);
            return 0;
        }

    if (maxStep >= 0.0)
        if (CVodeSetMaxStep(data->cvode_mem, maxStep) != CV_SUCCESS) {
            throw_java_exception(env, "jni_initialize: Failed to set maximum step size\n");
            // Free integrator memory
            CVodeFree(&data->cvode_mem);
            // Delete global java references
            (*env)->DeleteGlobalRef(env, data->ode);
            if (data->ef != NULL)
                (*env)->DeleteGlobalRef(env, data->ef);
            jni_free(env, data);
            return 0;
        }

    /* Specify user-data pointer for CVode */
    if (CVodeSetUserData(data->cvode_mem, (void*)data) != CV_SUCCESS) {
        throw_java_exception(env, "jni_initialize: Failed to set user-data\n");
        // Free integrator memory
        CVodeFree(&data->cvode_mem);
        // Delete global java references
        (*env)->DeleteGlobalRef(env, data->ode);
        if (data->ef != NULL)
            (*env)->DeleteGlobalRef(env, data->ef);
        jni_free(env, data);
        return 0;
    }

    /* Create serial vector of length neq */
    /*void* ydata = jni_malloc(env, sizeof(double) * neq);
    if (ydata == NULL) {
        throw_java_exception(env, "jni_initialize: Failed to allocate ydata with jni_malloc\n");
        // Free integrator memory
        CVodeFree(&data->cvode_mem);
        // Delete global java references
        (*env)->DeleteGlobalRef(env, data->ode);
        if (data->ef != NULL)
            (*env)->DeleteGlobalRef(env, data->ef);
        jni_free(env, data);
        return 0;
    }*/
    data->y = N_VMake_Serial(neq, data->x);
    if (data->y == NULL) {
        throw_java_exception(env, "jni_initialize: Failed to create serial vector y\n");
        //jni_free(env, ydata);
        // Free integrator memory
        CVodeFree(&data->cvode_mem);
        // Delete global java references
        (*env)->DeleteGlobalRef(env, data->ode);
        if (data->ef != NULL)
            (*env)->DeleteGlobalRef(env, data->ef);
        jni_free(env, data);
        return 0;
    }
    data->yTmp = N_VMake_Serial(neq, data->xTmp);
    if (data->yTmp == NULL) {
        throw_java_exception(env, "jni_initialize: Failed to create serial vector yTmp\n");
        //jni_free(env, ydata);
        // Free y vector
        N_VDestroy_Serial(data->y);
        // Free integrator memory
        CVodeFree(&data->cvode_mem);
        // Delete global java references
        (*env)->DeleteGlobalRef(env, data->ode);
        if (data->ef != NULL)
            (*env)->DeleteGlobalRef(env, data->ef);
        jni_free(env, data);
        return 0;
    }

    /* Call CVodeInit to initialize the integrator memory and specify the
    * user's right hand side function in y'=f(t,y), the inital time T0, and
    * the initial dependent variable vector y. */
    // For now use a default starting time t0
    if (CVodeInit(data->cvode_mem, &f_bridge, 0.0, data->y) != CV_SUCCESS) {
        throw_java_exception(env, "jni_initialize: Failed to initialize CVode problem\n");
        //jni_free(env, ydata);
        // Free y vector
        N_VDestroy_Serial(data->y);
        // Free yTmp vector
        N_VDestroy_Serial(data->yTmp);
        // Free integrator memory
        CVodeFree(&data->cvode_mem);
        // Delete global java references
        (*env)->DeleteGlobalRef(env, data->ode);
        if (data->ef != NULL)
            (*env)->DeleteGlobalRef(env, data->ef);
        jni_free(env, data);
        return 0;
    }

    /* Call CVodeSStolerances to specify the scalar relative tolerance
    * and scalar absolute tolerance */
    if (CVodeSStolerances(data->cvode_mem, relTol, absTol) != CV_SUCCESS) {
        throw_java_exception(env, "jni_initialize: Failed to set tolerances\n");
        //jni_free(env, ydata);
        // Free y vector
        N_VDestroy_Serial(data->y);
        // Free yTmp vector
        N_VDestroy_Serial(data->yTmp);
        // Free integrator memory
        CVodeFree(&data->cvode_mem);
        // Delete global java references
        (*env)->DeleteGlobalRef(env, data->ode);
        if (data->ef != NULL)
            (*env)->DeleteGlobalRef(env, data->ef);
        jni_free(env, data);
        return 0;
    }

    CALLBACK_G_PTR g_ptr;
    if (data->ng > 0)
        g_ptr = &g_bridge;
    else
        g_ptr = NULL;

    if (data->ng > 0) {
        /* Call CVodeRootInit to specify the root function g with data->ng components */
        if (CVodeRootInit(data->cvode_mem, data->ng, g_ptr) != CV_SUCCESS) {
            throw_java_exception(env, "jni_initialize: Failed to initialize root finding\n");
            //jni_free(env, ydata);
            // Free y vector
            N_VDestroy_Serial(data->y);
        // Free yTmp vector
        N_VDestroy_Serial(data->yTmp);
            // Free integrator memory
            CVodeFree(&data->cvode_mem);
            // Delete global java references
            (*env)->DeleteGlobalRef(env, data->ode);
            if (data->ef != NULL)
                (*env)->DeleteGlobalRef(env, data->ef);
            jni_free(env, data);
            return 0;
        }
    }

    /* Call CVDense to specify the CVDENSE dense linear solver */
    if (CVDense(data->cvode_mem, neq) != CV_SUCCESS) {
        throw_java_exception(env, "jni_initialize: Failed to initialize dense linear solver\n");
        //jni_free(env, ydata);
        // Free y vector
        N_VDestroy_Serial(data->y);
        // Free yTmp vector
        N_VDestroy_Serial(data->yTmp);
        // Free integrator memory
        CVodeFree(&data->cvode_mem);
        // Delete global java references
        (*env)->DeleteGlobalRef(env, data->ode);
        if (data->ef != NULL)
            (*env)->DeleteGlobalRef(env, data->ef);
        jni_free(env, data);
        return 0;
    }

    /* Specify user-defined Jacobian function for the dense linear solver */
    /*if (CVDlsSetDenseJacFn(data->cvode mem, jac_bridge) != CV_SUCCESS) {
        throw_java_exception(env, "jni_initialize: Failed to set Jacobian function\n");
        //jni_free(env, ydata);
        // Free y vector
        N_VDestroy_Serial(data->y);
        // Free yTmp vector
        N_VDestroy_Serial(data->yTmp);
        // Free integrator memory
        CVodeFree(&data->cvode_mem);
        // Delete global java references
        (*env)->DeleteGlobalRef(env, data->ode);
        if (data->ef != NULL)
            (*env)->DeleteGlobalRef(env, data->ef);
        jni_free(env, data);
        return 0;
    }*/

    /*if (CVodeSetMaxNumSteps(data->cvode_mem, 1000) != CV_SUCCESS) {
        throw_java_exception(env, "jni_initialize: Failed to set maximum number of steps\n");
        //jni_free(env, ydata);
        // Free y vector
        N_VDestroy_Serial(data->y);
        // Free yTmp vector
        N_VDestroy_Serial(data->yTmp);
        // Free integrator memory
        CVodeFree(&data->cvode_mem);
        // Delete global java references
        (*env)->DeleteGlobalRef(env, data->ode);
        if (data->ef != NULL)
            (*env)->DeleteGlobalRef(env, data->ef);
        jni_free(env, data);
        return 0;
    }*/

    if (data->ng > 0) {
        /* Create rootsfound array of length ng */
        data->rootsfound = jni_malloc(env, sizeof(int) * ng);
        if (data->rootsfound == NULL) {
            throw_java_exception(env, "jni_initialize: Failed to allocate rootsfound array with jni_malloc\n");
            //jni_free(env, ydata);
            // Free yTmp vector
            N_VDestroy_Serial(data->yTmp);
            // Free y vector
            N_VDestroy_Serial(data->y);
            // Free integrator memory
            CVodeFree(&data->cvode_mem);
            // Delete global java references
            (*env)->DeleteGlobalRef(env, data->ode);
            if (data->ef != NULL)
                (*env)->DeleteGlobalRef(env, data->ef);
            jni_free(env, data);
            return 0;
        }
    } else
        data->rootsfound = NULL;

    // Pass pointer-address of cvode_solver_data back to Java
    jlong jni_pointer = (jlong)data;
    return jni_pointer;
}

JNIEXPORT void JNICALL Java_ch_ethz_bhepp_ode_cvode_CVodeSolver_jni_1dispose
  (JNIEnv *env, jobject obj, jlong jni_pointer) {
    struct cvode_solver_data* data = (struct cvode_solver_data*)jni_pointer;
    if (data != NULL) {
        // Free rootsfound data
        if (data->rootsfound != NULL)
            jni_free(env, data->rootsfound);
        // Free yTmp vector
        N_VDestroy_Serial(data->yTmp);
        // Free y vector data
        //jni_free(env, NV_DATA_S(data->y));
        // Free y vector
        if (data->y != NULL)
            N_VDestroy_Serial(data->y);
        // Free integrator memory
        if (data->cvode_mem != NULL)
            CVodeFree(&data->cvode_mem);
        // Delete global java references
        if (data->ode != NULL)
            (*env)->DeleteGlobalRef(env, data->ode);
        if (data->ef != NULL)
            (*env)->DeleteGlobalRef(env, data->ef);
        // Free data structure
        jni_free(env, data);
    }
}

JNIEXPORT void JNICALL Java_ch_ethz_bhepp_ode_cvode_CVodeSolver_jni_1reinitialize
  (JNIEnv *env, jobject obj, jlong jni_pointer, jdouble t) {
    struct cvode_solver_data* data = (struct cvode_solver_data*)jni_pointer;
    if (CVodeReInit(data->cvode_mem, t, data->y) != CV_SUCCESS) {
        throw_java_exception(env, "jni_reinitialize: Failed to reinitialize CVode\n");
        return;
    }
}
/*JNIEXPORT void JNICALL Java_ch_ethz_bhepp_ode_cvode_CVodeSolver_jni_1quick_1reinitialize
  (JNIEnv *env, jobject obj, jlong jni_pointer, jdouble t) {
    struct cvode_solver_data* data = (struct cvode_solver_data*)jni_pointer;
    CVodeMem cv_mem = (CVodeMem)data->cvode_mem;
    cv_mem->cv_tn = t;
    N_VScale(1.0, data->y, cv_mem->cv_zn[0]);
}*/
JNIEXPORT void JNICALL Java_ch_ethz_bhepp_ode_cvode_CVodeSolver_jni_1reinitializeOneStep
  (JNIEnv *env, jobject obj, jlong jni_pointer, jdouble t) {
    struct cvode_solver_data* data = (struct cvode_solver_data*)jni_pointer;
    if (CVodeReInit(data->cvode_mem, t, data->y) != CV_SUCCESS) {
        throw_java_exception(env, "jni_reinitialize: Failed to reinitialize CVode\n");
        return;
    }
}

JNIEXPORT jdouble JNICALL Java_ch_ethz_bhepp_ode_cvode_CVodeSolver_jni_1integrate
  (JNIEnv *env, jobject obj, jlong jni_pointer, jdouble t1) {
    struct cvode_solver_data* data = (struct cvode_solver_data*)jni_pointer;
    int flag;
    // Reset event indicator
    data->eventIndex = -1;
    double tout;
    // Call CVode solver
    flag = CVode(data->cvode_mem, t1, data->y, &tout, CV_NORMAL);

    if (flag == CV_ROOT_RETURN) {
        // Get event index
        if (CVodeGetRootInfo(data->cvode_mem, data->rootsfound) != CV_SUCCESS) {
            throw_java_exception_error_code(env, "jni_integrate: Failed to get event index\n", flag);
            return NAN;
        }
        // Report one event to the solver
        int i;
        for (i=0; i < data->ng; i++)
            if (data->rootsfound[i] != 0) {
                data->eventIndex = i;
                break;
            }
    }

    if (flag < 0) {
        if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
            fprintf(stderr, "jni_integrate: Java exception was thrown during ODE integration\n");
            return NAN;
        }
        throw_java_exception_error_code(env, "jni_integrate: CVode integration failed\n", flag);
        return NAN;
    }

    return tout;
}

JNIEXPORT jdouble JNICALL Java_ch_ethz_bhepp_ode_cvode_CVodeSolver_jni_1integrateOneStep
  (JNIEnv *env, jobject obj, jlong jni_pointer, jdouble t1) {
    struct cvode_solver_data* data = (struct cvode_solver_data*)jni_pointer;
    int flag;
    // Reset event indicator
    data->eventIndex = -1;
    double tout;
    // Call CVode solver
    flag = CVode(data->cvode_mem, t1, data->y, &tout, CV_ONE_STEP);
    if (flag < 0) {
        if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
            fprintf(stderr, "jni_integrate: Java exception was thrown during ODE integration\n");
            return NAN;
        }
        throw_java_exception_error_code(env, "jni_integrate: CVode integration failed\n", flag);
        return NAN;
    }

    if (flag == CV_ROOT_RETURN) {
        // Get event index
        if (CVodeGetRootInfo(data->cvode_mem, data->rootsfound) != CV_SUCCESS) {
            throw_java_exception_error_code(env, "jni_integrate: Failed to get event index\n", flag);
            return NAN;
        }
        // Report one event to the solver
        int i;
        for (i=0; i < data->ng; i++)
            if (data->rootsfound[i] != 0) {
                data->eventIndex = i;
                break;
            }
    }

    return tout;
}

JNIEXPORT void JNICALL Java_ch_ethz_bhepp_ode_cvode_CVodeSolver_jni_1computeInterpolatedSolution
  (JNIEnv *env, jobject obj, jlong jni_pointer, jdouble t) {
    struct cvode_solver_data* data = (struct cvode_solver_data*)jni_pointer;
    int flag;
    // Compute interpolated solution
    int k = 0;
    flag = CVodeGetDky(data->cvode_mem, t, k, data->yTmp);
    if (flag < 0) {
        throw_java_exception_error_code(env, "jni_integrate: CVode interpolation failed\n", flag);
        return;
    }
}

JNIEXPORT jint JNICALL Java_ch_ethz_bhepp_ode_cvode_CVodeSolver_jni_1getEventIndex
  (JNIEnv *env, jobject obj, jlong jni_pointer) {
    struct cvode_solver_data* data = (struct cvode_solver_data*)jni_pointer;
    return data->eventIndex;
}

JNIEXPORT void JNICALL Java_ch_ethz_bhepp_ode_cvode_CVodeSolver_jni_1resetEventIndex
  (JNIEnv *env, jobject obj, jlong jni_pointer) {
    struct cvode_solver_data* data = (struct cvode_solver_data*)jni_pointer;
    data->eventIndex = -1;
}

/*
 * Class:     ch_ethz_bhepp_ode_cvode_CVodeSolver
 * Method:    jni_getCurrentStep
 * Signature: (J)D
 */
JNIEXPORT jdouble JNICALL Java_ch_ethz_bhepp_ode_cvode_CVodeSolver_jni_1getCurrentStep
  (JNIEnv *env, jobject obj, jlong jni_pointer) {
    struct cvode_solver_data* data = (struct cvode_solver_data*)jni_pointer;
    jdouble hcur;
    int flag;
    flag = CVodeGetCurrentStep(data->cvode_mem, &hcur);
    if (flag < 0) {
        throw_java_exception_error_code(env, "jni_integrate: CVodeGetCurrentStep failed\n", flag);
        return NAN;
    }
    return hcur;
}

/*
 * Class:     ch_ethz_bhepp_ode_cvode_CVodeSolver
 * Method:    jni_getLastStep
 * Signature: (J)D
 */
JNIEXPORT jdouble JNICALL Java_ch_ethz_bhepp_ode_cvode_CVodeSolver_jni_1getLastStep
  (JNIEnv *env, jobject obj, jlong jni_pointer) {
    struct cvode_solver_data* data = (struct cvode_solver_data*)jni_pointer;
    jdouble hlast;
    int flag;
    flag = CVodeGetLastStep(data->cvode_mem, &hlast);
    if (flag < 0) {
        throw_java_exception_error_code(env, "jni_integrate: CVodeGetLastStep failed\n", flag);
        return NAN;
    }
    return hlast;
}

/*
 * Class:     ch_ethz_bhepp_ode_cvode_CVodeSolver
 * Method:    jni_setInitStep
 * Signature: (JD)V
 */
JNIEXPORT void JNICALL Java_ch_ethz_bhepp_ode_cvode_CVodeSolver_jni_1setInitStep
  (JNIEnv *env, jobject obj, jlong jni_pointer, jdouble hin) {
    struct cvode_solver_data* data = (struct cvode_solver_data*)jni_pointer;
    int flag;
    flag = CVodeSetInitStep(data->cvode_mem, hin);
    if (flag < 0) {
        throw_java_exception_error_code(env, "jni_integrate: CVodeSetInitStep failed\n", flag);
        return;
    }
    return;
}

/*
 * Class:     ch_ethz_bhepp_ode_cvode_CVodeSolver
 * Method:    jni_setMinStep
 * Signature: (JD)V
 */
JNIEXPORT void JNICALL Java_ch_ethz_bhepp_ode_cvode_CVodeSolver_jni_1setMinStep
  (JNIEnv *env, jobject obj, jlong jni_pointer, jdouble hmin) {
    struct cvode_solver_data* data = (struct cvode_solver_data*)jni_pointer;
    int flag;
    flag = CVodeSetMinStep(data->cvode_mem, hmin);
    if (flag < 0) {
        throw_java_exception_error_code(env, "jni_integrate: CVodeSetMinStep failed\n", flag);
        return;
    }
    return;
}

/*
 * Class:     ch_ethz_bhepp_ode_cvode_CVodeSolver
 * Method:    jni_setMaxStep
 * Signature: (JD)V
 */
JNIEXPORT void JNICALL Java_ch_ethz_bhepp_ode_cvode_CVodeSolver_jni_1setMaxStep
  (JNIEnv *env, jobject obj, jlong jni_pointer, jdouble hmax) {
    struct cvode_solver_data* data = (struct cvode_solver_data*)jni_pointer;
    int flag;
    flag = CVodeSetMaxStep(data->cvode_mem, hmax);
    if (flag < 0) {
        throw_java_exception_error_code(env, "jni_integrate: CVodeSetMaxStep failed\n", flag);
        return;
    }
    return;
}
