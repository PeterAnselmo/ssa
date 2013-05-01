#ifndef UTILS_CPP
#define UTILS_CPP

void check_ptr(char* ptr){
    if(ptr == NULL){
        printf("Memory Allocation Error.\n");
        exit(1);
    }
}

#endif
