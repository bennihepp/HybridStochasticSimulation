#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/dense_output_runge_kutta.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include "ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver.h"

#include "ch_ethz_bhepp_ode_boost_NativeStepperType.h"
#include "Utilities.hpp"

using namespace boost::numeric::odeint;

typedef boost::numeric::ublas::vector<double> vector_type;
typedef double time_type;
typedef vector_type state_type;


template< class ode_system >
class stepper_wrapper {
public:

	virtual ~stepper_wrapper() { };

	virtual void initialize(const state_type& x0, const time_type t0, const time_type dt0) = 0;

	virtual std::pair< time_type, time_type > do_step(ode_system& sys) = 0;

	virtual void calc_state(const time_type t_inter, state_type& x) = 0;

	virtual const time_type current_time() = 0;

	virtual const state_type& current_state() = 0;

	virtual const time_type current_time_step() = 0;

};

template < typename base_stepper_type, typename ode_system >
class stepper_wrapper_impl : public stepper_wrapper< ode_system > {
	base_stepper_type base_stepper;

public:

	stepper_wrapper_impl(base_stepper_type base_stepper)
	: base_stepper(base_stepper) { }

	void initialize(const state_type& x0, const time_type t0, const time_type dt0) {
		base_stepper.initialize(x0, t0, dt0);
	}

	std::pair< time_type, time_type > do_step(ode_system& sys) {
		return base_stepper.do_step(sys);
	}

	void calc_state(const time_type t_inter, state_type& x) {
		base_stepper.calc_state(t_inter, x);
	}

	const time_type current_time() {
		return base_stepper.current_time();
	}

	const state_type& current_state() {
		return base_stepper.current_state();
	}

	const time_type current_time_step() {
		return base_stepper.current_time_step();
	}

};

struct ode_system_wrapper;

struct boost_solver_data {

//	dense_output_runge_kutta< runge_kutta_dopri5< state_type > > stepper;
//	runge_kutta_dopri5< state_type > stepper;
	stepper_wrapper< ode_system_wrapper >* stepper;
	ode_system_wrapper* ode_system;
//	dense_dopri5_stepper_type stepper;
    JNIEnv* env;
    jobject ode;
    jmethodID vectorFieldMethodID;
    jdouble initialStepSize;
    jdouble t;
    jdouble t1;
    jdouble* x;
    jdouble* xTmp;
    jdouble* xDot;
    state_type y;
    state_type yTmp;
    state_type yDot;
	jdouble lastStepSize;
	jdouble currentStepSize;

	void operator()(const state_type &x, state_type &dxdt, double t) {
//	    for (state_type::const_iterator it=x.begin(); it != x.end(); ++it) {
		for (int i=0; i < x.size(); i++)
	    	xTmp[i] = x[i];
	    env->CallVoidMethod(ode, vectorFieldMethodID, t);
	    if (env->ExceptionCheck() == JNI_TRUE) {
	        fprintf(stderr, "f_direct_bridge: Failed to call computeVectorField!\n");
	        return;
	    }
//	    for (state_type::const_iterator it=x.begin(); it != x.end(); ++it)
		for (int i=0; i < dxdt.size(); i++)
	    	dxdt[i] = xDot[i];
    }

};

struct ode_system_wrapper {
	boost_solver_data* data;

	ode_system_wrapper(boost_solver_data* data)
	: data(data) { }

	void operator()(const state_type &x, state_type &dxdt, double t) {
		data->operator()(x, dxdt, t);
	}
};

typedef result_of::make_dense_output<
    runge_kutta_dopri5< state_type > >::type dense_dopri5_stepper_type;
typedef bulirsch_stoer_dense_out< state_type > dense_bulirsch_stoer_stepper_type;
//typedef result_of::make_dense_output<
//    bulirsch_stoer< state_type > >::type dense_bulirsch_stoer_stepper_type;
//typedef result_of::make_dense_output<
//	rosenbrock4< double > >::type dense_rosenbrock4_stepper_type;

template class stepper_wrapper_impl< dense_dopri5_stepper_type, ode_system_wrapper >;
template class stepper_wrapper_impl< dense_bulirsch_stoer_stepper_type, ode_system_wrapper >;
//template class stepper_wrapper_impl< dense_rosenbrock4_stepper_type, ode_system_wrapper >;

/*
 * Class:     ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver
 * Method:    jni_initialize
 * Signature: (Lch/ethz/bhepp/ode/BufferOdeAdapter;Ljava/nio/DoubleBuffer;Ljava/nio/DoubleBuffer;Ljava/nio/DoubleBuffer;DDDI)J
 */
JNIEXPORT jlong JNICALL Java_ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver_jni_1initialize
	(JNIEnv *env, jobject obj, jobject ode, jobject xBuffer, jobject xTmpBuffer, jobject xDotBuffer,
			jdouble initialStepSize, jdouble relTol, jdouble absTol, jint stepperType) {
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

    switch (stepperType) {
    case ch_ethz_bhepp_ode_boost_NativeStepperType_BOOSTSTEPPERTYPE_DORMANDPRINCE5:
    {
    	dense_dopri5_stepper_type base_stepper1 = make_dense_output< runge_kutta_dopri5< state_type > >(absTol, relTol);
    	data->stepper = new stepper_wrapper_impl< dense_dopri5_stepper_type, ode_system_wrapper >(base_stepper1);
    	break;
    }
    case ch_ethz_bhepp_ode_boost_NativeStepperType_BOOSTSTEPPERTYPE_BULIRSCHSTOER:
    {
    	dense_bulirsch_stoer_stepper_type base_stepper2 = dense_bulirsch_stoer_stepper_type(absTol, relTol);
    	data->stepper = new stepper_wrapper_impl< dense_bulirsch_stoer_stepper_type, ode_system_wrapper >(base_stepper2);
    	break;
    }
//    case ch_ethz_bhepp_ode_boost_NativeStepperType_BOOSTSTEPPERTYPE_ROSENBROCK4:
//    {
//    	dense_rosenbrock4_stepper_type base_stepper3 = make_dense_output< rosenbrock4< double > >(absTol, relTol);
//    	data->stepper = new stepper_wrapper_impl< dense_rosenbrock4_stepper_type, ode_system_wrapper >(base_stepper3);
//    	break;
//    }
    default:
	{
        throw_java_exception(env, "Unknown stepper type\n");
    	delete data;
    	return 0;
	}
    }

    // FIXME
    data->ode_system = new ode_system_wrapper(data);
    data->initialStepSize = initialStepSize;

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
    data->yTmp.resize(dimension, false);
    data->yDot.resize(dimension, false);

    data->lastStepSize = data->initialStepSize;
    data->currentStepSize = data->initialStepSize;

    // Pass pointer-address of boost_solver_data back to Java
    jlong jni_pointer = (jlong)data;
    return jni_pointer;
}

/*
 * Class:     ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver
 * Method:    jni_dispose
 * Signature: (J)V
 */
JNIEXPORT void JNICALL Java_ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver_jni_1dispose
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

/*
 * Class:     ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver
 * Method:    jni_prepareStep
 * Signature: (JDD)V
 */
JNIEXPORT void JNICALL Java_ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver_jni_1prepareStep
	(JNIEnv *env, jobject obj, jlong jni_pointer, jdouble t0, jdouble t1) {
    boost_solver_data* data = reinterpret_cast< boost_solver_data* >(jni_pointer);
    data->t = t0;
    data->t1 = t1;
    for (int i=0; i < data->y.size(); i++)
    	data->y[i] = data->x[i];
    data->stepper->initialize(data->y, data->t, data->currentStepSize);
}

/*
 * Class:     ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver
 * Method:    jni_getCurrentState
 * Signature: (J)D
 */
JNIEXPORT jdouble JNICALL Java_ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver_jni_1getCurrentState
  (JNIEnv *env, jobject obj, jlong jni_pointer) {
    boost_solver_data* data = reinterpret_cast< boost_solver_data* >(jni_pointer);
    const dense_dopri5_stepper_type::state_type& current_state = data->stepper->current_state();
    for (int i=0; i < data->y.size(); i++)
    	data->xTmp[i] = current_state[i];
    return data->stepper->current_time();
}

/*
 * Class:     ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver
 * Method:    jni_integrateStep
 * Signature: (J)D
 */
JNIEXPORT jdouble JNICALL Java_ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver_jni_1integrateStep
(JNIEnv *env, jobject obj, jlong jni_pointer) {
    boost_solver_data* data = reinterpret_cast< boost_solver_data* >(jni_pointer);
    std::pair< double , double > times = data->stepper->do_step(*data->ode_system);
    data->lastStepSize = data->currentStepSize;
    data->currentStepSize = data->stepper->current_time_step();
    return times.second;
}

/*
 * Class:     ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver
 * Method:    jni_computeInterpolatedSolution
 * Signature: (JD)V
 */
JNIEXPORT void JNICALL Java_ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver_jni_1computeInterpolatedSolution
  (JNIEnv *env, jobject obj, jlong jni_pointer, jdouble t) {
    boost_solver_data* data = reinterpret_cast< boost_solver_data* >(jni_pointer);
    data->stepper->calc_state(t, data->yTmp);
    for (int i=0; i < data->yTmp.size(); i++)
    	data->xTmp[i] = data->yTmp[i];
}

/*
 * Class:     ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver
 * Method:    jni_setCurrentStepSize
 * Signature: (JD)V
 */
JNIEXPORT void JNICALL Java_ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver_jni_1setCurrentStepSize
  (JNIEnv *env, jobject obj, jlong jni_pointer, jdouble stepSize) {
    boost_solver_data* data = reinterpret_cast< boost_solver_data* >(jni_pointer);
    data->currentStepSize = stepSize;
}

/*
 * Class:     ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver
 * Method:    jni_getCurrentStepSize
 * Signature: (J)D
 */
JNIEXPORT jdouble JNICALL Java_ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver_jni_1getCurrentStepSize
  (JNIEnv *env, jobject obj, jlong jni_pointer) {
    boost_solver_data* data = reinterpret_cast< boost_solver_data* >(jni_pointer);
    return data->currentStepSize;
}

/*
 * Class:     ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver
 * Method:    jni_getLastStepSize
 * Signature: (J)D
 */
JNIEXPORT jdouble JNICALL Java_ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver_jni_1getLastStepSize
(JNIEnv *env, jobject obj, jlong jni_pointer) {
  boost_solver_data* data = reinterpret_cast< boost_solver_data* >(jni_pointer);
  return data->lastStepSize;
}
