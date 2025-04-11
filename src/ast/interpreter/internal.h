# 2 "ast/interpreter/internal.h"

/**
# Internally defined functions */

void * malloc (long size) {}
void * calloc (long nmemb, long size) {}
void * realloc (void * ptr, long size) {}
void free (void * ptr) {}
void * memset (void * s, int c, long n) {}
void * memcpy (void * dest, const void * src, long n) {}
char * strdup (const char * s) {}
char * strcpy (char *dest, const char *src) {}
char * strcat (char * dest, const char * src) {}
long strlen (const char * s) {}
int strcmp (const char * s1, const char * s2) {}
int strncmp (const char * s1, const char * s2, long size) {}

double fabs (double x) {}
double sqrt (double x) {}
double exp (double a) {}
double log (double x) {}
double log10 (double a) {}
double sin (double x) {}
double cos (double a) {}
double tan (double x) {}
double asin (double x) {}
double acos (double a) {}
double atan (double x) {}
double sinh (double a) {}
double cosh (double a) {}
double tanh (double a) {}
double asinh (double a) {}
double acosh (double a) {}
double atanh (double a) {}

double atan2 (double y, double x) {}
double pow (double x, double y) {}

void interpreter_verbosity (int verbosity) {}
void interpreter_maximum_iterations (int maximum_iterations) {}
void display_value (void * value) {}
void reset_field_value (real * field, const char * name, real val, int block) {}
