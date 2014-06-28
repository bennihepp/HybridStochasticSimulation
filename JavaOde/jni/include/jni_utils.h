#ifndef _INCLUDED_JNI_UTILS_H_
#define _INCLUDED_JNI_UTILS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <jni.h>

#include <memory.h>

extern void* jni_malloc(JNIEnv* env, size_t size);

extern void jni_free(JNIEnv* env, void* ptr);

#ifdef __cplusplus
}
#endif

#endif

