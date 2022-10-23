#ifndef TESTDRIVE_ERRORS_H
#define TESTDRIVE_ERRORS_H

#define ERR_INVALID_ARG_NUM 1
#define ERR_INVALID_IN_FILE 2
#define ERR_MEM_ALLOC_FAILED 3
#define ERR_READ_WRITE_FAILED 4
#define ERR_FSEEK_FAILED 5
#define ERR_DIVISION_BY_ZERO 6

void error_msg(int err_code);

#endif
