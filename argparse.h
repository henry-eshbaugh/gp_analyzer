#ifndef _argparse_h
#define _argparse_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ARGV_COUNTER _argv_counter
#define PREV_ARG ARGV_COUNTER--
#define NEXT_ARG ARGV_COUNTER++
#define THIS_ARG argv[ARGV_COUNTER]
#define NEXT_ARG_MUST_EXIST do { NEXT_ARG; CHECK_ARGV; } while (0)
#define LOOP_ARGV for (int ARGV_COUNTER = 1; ARGV_COUNTER < argc; NEXT_ARG)
#define IF_ARG(x) if (!strcasecmp(argv[ARGV_COUNTER], x))
#define IF_NOT_ARG(x) if (strcasecmp(argv[ARGV_COUNTER], x))
#define IF_ARG_EXISTS if (!(ARGV_COUNTER < argc)) else
#define IF_ARG_NOT_EXISTS if (!(ARGV_COUNTER >= argc)); else usage()
#define CHECK_ARGV if (!(ARGV_COUNTER >= argc)); else usage()

#define CHECK_OK(cond, ...) if (cond) { fprintf(stderr, __VA_ARGS__); usage(); }

#endif /* !defined _argparse_h */
