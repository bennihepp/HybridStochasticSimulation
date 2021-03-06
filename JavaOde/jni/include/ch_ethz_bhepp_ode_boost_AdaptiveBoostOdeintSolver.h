/* DO NOT EDIT THIS FILE - it is machine generated */
#include <jni.h>
/* Header for class ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver */

#ifndef _Included_ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver
#define _Included_ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver
#ifdef __cplusplus
extern "C" {
#endif
/*
 * Class:     ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver
 * Method:    jni_initialize
 * Signature: (Lch/ethz/bhepp/ode/BufferOdeAdapter;Ljava/nio/DoubleBuffer;Ljava/nio/DoubleBuffer;Ljava/nio/DoubleBuffer;DDDI)J
 */
JNIEXPORT jlong JNICALL Java_ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver_jni_1initialize
  (JNIEnv *, jobject, jobject, jobject, jobject, jobject, jdouble, jdouble, jdouble, jint);

/*
 * Class:     ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver
 * Method:    jni_dispose
 * Signature: (J)V
 */
JNIEXPORT void JNICALL Java_ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver_jni_1dispose
  (JNIEnv *, jobject, jlong);

/*
 * Class:     ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver
 * Method:    jni_prepareStep
 * Signature: (JDD)V
 */
JNIEXPORT void JNICALL Java_ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver_jni_1prepareStep
  (JNIEnv *, jobject, jlong, jdouble, jdouble);

/*
 * Class:     ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver
 * Method:    jni_setCurrentStepSize
 * Signature: (JD)V
 */
JNIEXPORT void JNICALL Java_ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver_jni_1setCurrentStepSize
  (JNIEnv *, jobject, jlong, jdouble);

/*
 * Class:     ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver
 * Method:    jni_getCurrentState
 * Signature: (J)D
 */
JNIEXPORT jdouble JNICALL Java_ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver_jni_1getCurrentState
  (JNIEnv *, jobject, jlong);

/*
 * Class:     ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver
 * Method:    jni_integrateStep
 * Signature: (J)D
 */
JNIEXPORT jdouble JNICALL Java_ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver_jni_1integrateStep
  (JNIEnv *, jobject, jlong);

/*
 * Class:     ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver
 * Method:    jni_computeInterpolatedSolution
 * Signature: (JD)V
 */
JNIEXPORT void JNICALL Java_ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver_jni_1computeInterpolatedSolution
  (JNIEnv *, jobject, jlong, jdouble);

/*
 * Class:     ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver
 * Method:    jni_getCurrentStepSize
 * Signature: (J)D
 */
JNIEXPORT jdouble JNICALL Java_ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver_jni_1getCurrentStepSize
  (JNIEnv *, jobject, jlong);

/*
 * Class:     ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver
 * Method:    jni_getLastStepSize
 * Signature: (J)D
 */
JNIEXPORT jdouble JNICALL Java_ch_ethz_bhepp_ode_boost_AdaptiveBoostOdeintSolver_jni_1getLastStepSize
  (JNIEnv *, jobject, jlong);

#ifdef __cplusplus
}
#endif
#endif
