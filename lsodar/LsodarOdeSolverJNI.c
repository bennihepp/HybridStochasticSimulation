#include <stdio.h>

#include "LsodarOdeSolver.h"

JNIEXPORT jint JNICALL Java_LsodarOdeSolver_callbackTest
  (JNIEnv* env, jclass class, jdoubleArray ja) {

    jsize n = (*env)->GetArrayLength(env, ja);
    jdouble *a = (*env)->GetDoubleArrayElements(env, ja, 0);

    jmethodID methodID = (*env)->GetStaticMethodID(env, class, "callback", "(D)D");

    a[n-1] = (*env)->CallStaticDoubleMethod(env, class, methodID, a[n-1]);

    (*env)->ReleaseDoubleArrayElements(env, ja, a, 0);

    return 0;
}
