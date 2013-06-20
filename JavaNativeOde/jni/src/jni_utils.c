#include "jni_utils.h"

struct Jalloc {
  jbyteArray jba;
  //jobject ref;
};

void* jni_malloc(JNIEnv* env, size_t size) {
    // Create a Java array that is big enough to hold a struct Jalloc
    // (for keeping a global object reference to the array) and the
    // requested amount of memory
    jbyteArray jba = (*env)->NewByteArray(env, (int)size + sizeof(struct Jalloc));
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "jni_malloc: Failed to create jba array\n");
        return NULL;
    }
    void* jbuffer = (void*)((*env)->GetByteArrayElements(env, jba, 0));
    if (jbuffer == NULL) {
        fprintf(stderr, "jni_malloc: Failed to get array elements of jba\n");
        return NULL;
    }
    struct Jalloc* pJalloc = (struct Jalloc*)jbuffer;
    //pJalloc->jba = jba;
    // Assign a global reference so byte array will persist until released
    pJalloc->jba = (*env)->NewGlobalRef(env, jba);
    if ((*env)->ExceptionCheck(env) == JNI_TRUE) {
        fprintf(stderr, "jni_malloc: Failed to create global reference to jba object\n");
        return NULL;
    }
    // Return address of the working area in jba
    return (void*)( ((char*)jbuffer) + sizeof(struct Jalloc));
}

void jni_free(JNIEnv* env, void* ptr) {
    if (ptr != NULL) {
        // Acquire the global reference of the byte array and release it
        void* buffer = (void*)( ((char*)ptr) - sizeof(struct Jalloc));
        struct Jalloc* pJalloc = (struct Jalloc*)buffer;
        jbyteArray jba = pJalloc->jba;
        (*env)->ReleaseByteArrayElements(env, pJalloc->jba, (jbyte*)buffer, JNI_ABORT);
        (*env)->DeleteGlobalRef(env, jba);
    }
}
