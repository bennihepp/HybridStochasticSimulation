#ifndef _INCLUDED_JNI_BOOST_UTILITIES_H_
#define _INCLUDED_JNI_BOOST_UTILITIES_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <jni.h>

#ifdef __cplusplus
}
#endif

int throw_java_exception(JNIEnv* env, const char* msg);

#endif
