#ifndef CUDASTRING_CU
#define CUDASTRING_CU

__device__ int cudaStrlen(const char* string){
    int length = 0;
    while(string[length++] != '\0');

    return length-1;
}

__device__ char* cudaStrncpy(char* dest, const char* source, int length){
    int pos = 0;
    while(pos < length && source[pos] != '\0'){
        dest[pos] = source[pos];
        ++pos;
    }
    while(pos < length){
        dest[pos] = '\0';
        ++pos;
    }
    return dest;
}
    


#endif
