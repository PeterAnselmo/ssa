#ifndef UTILS_CPP
#define UTILS_CPP

#define check_ptr(ptr) { _check_ptr((ptr), __FILE__, __LINE__); }
inline void _check_ptr(char* ptr, const char* file, int line){
    if(ptr == NULL){
        printf("Memory Allocation Error at %s:%d\n", file, line);
        exit(1);
    }
}

#endif
