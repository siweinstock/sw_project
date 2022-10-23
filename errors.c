#include <stdio.h>
#include <stdlib.h>
#include "errors.h"

void error_msg(int err_code) {
    char* msg;

    printf("ERROR: ");
    switch(err_code) {
        case ERR_INVALID_ARG_NUM:
            msg = "Invalid number of arguments.";
            break;
        case ERR_INVALID_IN_FILE:
            msg = "Input file does not exist or cannot be opened.";
            break;

        default:
            msg = "Unknown error.\n";
    }

    printf("%s\n", msg);
    exit(err_code);
}
