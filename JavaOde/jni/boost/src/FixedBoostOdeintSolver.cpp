#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/dense_output_runge_kutta.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include "ch_ethz_bhepp_ode_boost_FixedBoostOdeintSolver.h"

#include "Utilities.hpp"

using namespace boost::numeric::odeint;

typedef boost::numeric::ublas::vector<double> vector_type;
typedef vector_type state_type;
typedef result_of::make_dense_output<
    runge_kutta_dopri5< state_type > >::type dense_dopri5_stepper_type;

//typedef int (*CALLBACK_F_PTR)(double t, N_Vector y, N_Vector ydot, void* user_data);
//typedef int (*CALLBACK_G_PTR)(double t, N_Vector y, double* gout, void* user_data);


struct ode_system;

struct boost_solver_data {
//	dense_output_runge_kutta< runge_kutta_dopri5< state_type > > stepper;
//	runge_kutta_dopri5< state_type > stepper;
	ode_system* system;
	dense_dopri5_stepper_type stepper;
    JNIEnv* env;
    jobject ode;
    jmethodID vectorFieldMethodID;
    jdouble stepSize;
    jdouble t;
    jdouble t1;
    jdouble* x;
    jdouble* xTmp;
    jdouble* xDot;
    state_type y;
    state_type yDot;
};

struct ode_system {

    boost_solver_data* data;

	void operator()(const state_type &x, state_type &dxdt, double t )
    {
	    JNIEnv *env = data->env;
//	    for (state_type::const_iterator it=x.begin(); it != x.end(); ++it) {
		for (int i=0; i < x.size(); i++)
	    	data->xTmp[i] = x[i];
	    env->CallVoidMethod(data->ode, data->vectorFieldMethodID, t);
	    if (env->ExceptionCheck() == JNI_TRUE) {
	        fprintf(stderr, "f_direct_bridge: Failed to call computeVectorField!\n");
	        return;
	    }
//	    for (state_type::const_iterator it=x.begin(); it != x.end(); ++it)
		for (int i=0; i < dxdt.size(); i++)
	    	dxdt[i] = data->xDot[i];
    }

};

struct ode_system_wrapper {

	ode_system* system;

	ode_system_wrapper(ode_system* system)
		: system(system) { }
	void operator()(const state_type &x, state_type &dxdt, double t ) {
		system->operator ()(x, dxdt, t);
	}

};

/*int jac_bridge(long int N, double t, N Vector y, N Vector fy, DlsMat Jac,
               void *user data, N Vector tmp1, N Vector tmp2, N Vector tmp3) {
    return 0;
}*/

JNIEXPORT jlong JNICALL Java_ch_ethz_bhepp_ode_boost_FixedBoostOdeintSolver_jni_1initialize
	(JNIEnv *env, jobject obj, jobject ode, jobject xBuffer, jobject xTmpBuffer, jobject xDotBuffer,
			jdouble stepSize, jdouble relTol, jdouble absTol, jint stepperType) {
    jclass odeCls = env->GetObjectClass(ode);
    if (env->ExceptionCheck() == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to get class of ode object\n");
        return 0;
    }
    jmethodID methodID = env->GetMethodID(odeCls, "getDimensionOfVectorField", "()I");
    if (env->ExceptionCheck() == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to get method ID of getDimensionOfVectorField\n");
        return 0;
    }
    jint jdimension = env->CallIntMethod(ode, methodID);
    state_type::size_type dimension = jdimension;
    if (env->ExceptionCheck() == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to call getDimensionOfVectorField\n");
        return 0;
    }

    // Allocate and initialize boost_solver_data
    boost_solver_data* data = new boost_solver_data();
    // FIXME
    if (env->ExceptionCheck() == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to allocate boost_solver_data\n");
        return 0;
    }
    data->stepSize = stepSize;

    data->env = env;
    // Acquire method IDs of callbacks and put them into boost_solver_data
    data->ode = env->NewGlobalRef(ode);
    if (env->ExceptionCheck() == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to create global reference to ode object\n");
        delete data;
        return 0;
    }
    data->vectorFieldMethodID = env->GetMethodID(odeCls, "computeVectorField", "(D)V");
    if (env->ExceptionCheck() == JNI_TRUE) {
        fprintf(stderr, "jni_initialize: Failed to get method ID of computeVectorField\n");
        // Delete global java references
        env->DeleteGlobalRef(data->ode);
        delete data;
        return 0;
    }

    // Check capacities of direct buffer objects
    if (env->GetDirectBufferCapacity(xBuffer) < dimension) {
        throw_java_exception(env, "jni_initialize: Direct buffer xBuffer is not big enough\n");
        // Delete global java references
        env->DeleteGlobalRef(data->ode);
        delete data;
        return 0;
    }
    if (env->GetDirectBufferCapacity(xTmpBuffer) < dimension) {
        throw_java_exception(env, "jni_initialize: Direct buffer xTmpBuffer is not big enough\n");
        // Delete global java references
        env->DeleteGlobalRef(data->ode);
        delete data;
        return 0;
    }
    if (env->GetDirectBufferCapacity(xDotBuffer) < dimension) {
        throw_java_exception(env, "jni_initialize: Direct buffer xDotBuffer is not big enough\n");
        // Delete global java references
        env->DeleteGlobalRef(data->ode);
        delete data;
        return 0;
    }
    // Get addresses of direct buffer objects
    data->x = static_cast<jdouble*>(env->GetDirectBufferAddress(xBuffer));
    if (data->x == NULL) {
        fprintf(stderr, "jni_initialize: Failed to get direct address of xBuffer\n");
        // Delete global java references
        env->DeleteGlobalRef(data->ode);
        delete data;
        return 0;
    }
    data->xTmp = static_cast<jdouble*>(env->GetDirectBufferAddress(xTmpBuffer));
    if (data->xTmp == NULL) {
        fprintf(stderr, "jni_initialize: Failed to get direct address of xTmpBuffer\n");
        // Delete global java references
        env->DeleteGlobalRef(data->ode);
        delete data;
        return 0;
    }
    data->xDot = static_cast<jdouble*>(env->GetDirectBufferAddress(xDotBuffer));
    if (data->xDot == NULL) {
        fprintf(stderr, "jni_initialize: Failed to get direct address of xDotBuffer\n");
        // Delete global java references
        env->DeleteGlobalRef(data->ode);
        delete data;
        return 0;
    }

    data->y.resize(dimension, false);
    data->yDot.resize(dimension, false);

    data->system = new ode_system;
    data->system->data = data;

    // Pass pointer-address of boost_solver_data back to Java
    jlong jni_pointer = (jlong)data;
    return jni_pointer;
}

JNIEXPORT void JNICALL Java_ch_ethz_bhepp_ode_boost_FixedBoostOdeintSolver_jni_1dispose
  (JNIEnv *env, jobject obj, jlong jni_pointer) {
    boost_solver_data* data = (boost_solver_data*)jni_pointer;
    if (data != NULL) {
        // Delete global java references
        if (data->ode != NULL)
            env->DeleteGlobalRef(data->ode);
        // Free data structure
        delete data;
    }
}

JNIEXPORT void JNICALL Java_ch_ethz_bhepp_ode_boost_FixedBoostOdeintSolver_jni_1prepareStep
	(JNIEnv *env, jobject obj, jlong jni_pointer, jdouble t0, jdouble t1) {
    boost_solver_data* data = reinterpret_cast< boost_solver_data* >(jni_pointer);
    data->t = t0;
    data->t1 = t1;
    for (int i=0; i < data->y.size(); i++)
    	data->y[i] = data->x[i];
    data->stepper.initialize(data->y, data->t, data->stepSize);
}

/*
 * Class:     ch_ethz_bhepp_ode_boost_FixedBoostOdeintSolver
 * Method:    jni_integrateOneStep
 * Signature: (J)D
 */
JNIEXPORT jdouble JNICALL Java_ch_ethz_bhepp_ode_boost_FixedBoostOdeintSolver_jni_1integrateOneStep
(JNIEnv *env, jobject obj, jlong jni_pointer) {
    boost_solver_data* data = reinterpret_cast< boost_solver_data* >(jni_pointer);
    std::pair< double , double > times = data->stepper.do_step(ode_system_wrapper(data->system));
    const dense_dopri5_stepper_type::state_type& current_state = data->stepper.current_state();
    for (int i=0; i < data->y.size(); i++)
    	data->x[i] = current_state[i];
    return times.second;
}
