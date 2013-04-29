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

__device__ int cudaStrcmp(const char* s1, const char* s2){
    unsigned char uc1, uc2;
    /* Move s1 and s2 to the first differing characters 
    in each string, or the ends of the strings if they
    are identical.  */
    while (*s1 != '\0' && *s1 == *s2) {
        s1++;
        s2++;
    }
    /* Compare the characters as unsigned char and
    return the difference.  */
    uc1 = (*(unsigned char *) s1);
    uc2 = (*(unsigned char *) s2);
    return ((uc1 < uc2) ? -1 : (uc1 > uc2));
}
    


#endif
