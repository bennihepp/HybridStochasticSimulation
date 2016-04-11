#include "Utilities.hpp"

int throw_java_exception(JNIEnv* env, const char* msg) {
    jclass excpClass = env->FindClass("ch/ethz/bhepp/ode/boost/JniException");
    if (env->ExceptionCheck() == JNI_TRUE) {
        fprintf(stderr, "throw_java_exception: Failed to find class ch.ethz.bhepp.ode.boost.JniException\n");
        fprintf(stderr, "Original message: %s\n", msg);
        return -1;
    }
    jmethodID constructorMethodID = env->GetMethodID(excpClass, "<init>", "(Ljava/lang/String;)V");
    if (env->ExceptionCheck() == JNI_TRUE) {
        fprintf(stderr, "throw_java_exception: Failed to get method ID of constructor for ch.ethz.bhepp.ode.cvode.CVodeSolver.JniException\n");
        fprintf(stderr, "Original message: %s\n", msg);
        return -1;
    }
    jthrowable excp = reinterpret_cast<jthrowable>(env->NewObject(excpClass, constructorMethodID, msg));
    if (env->ExceptionCheck() == JNI_TRUE) {
        fprintf(stderr, "throw_java_exception: Failed to create ch.ethz.bhepp.ode.boost.JniException object\n");
        fprintf(stderr, "Original message: %s\n", msg);
        return -1;
    }
//    int result = env->ThrowNew(excpClass, msg);
    int result = env->Throw(excp);
    if (result != 0) {
        fprintf(stderr, "throw_java_exception: Failed to throw exception\n");
        fprintf(stderr, "Original message: %s\n", msg);
        return -1;
    }
    return 0;
}
