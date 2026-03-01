/**
# Parameter Management Library

This library provides utilities for reading and managing input parameters from
configuration files, creating output directories, and backing up files. It's
designed to work with Basilisk simulations and supports various parameter types
including integers, doubles, and arrays.


The configuration file should contain one parameter per line in the format:

~~~bash
key = value
~~~

'value' can be an int, double or array. Comments are allowed with # character

Example:

~~~bash
# parameter list

iterations = 1000
timestep = 0.01   # with comments
coefficients = [1.0,2.5,3.7]
~~~

## Global Variables and Data Structures
*/

char dpath[80]; // name of output dir
char file_param[80] = "params.in"; // name param file


// for mkdir
#include <sys/stat.h>
#include <sys/types.h>


/**
## String Utilities

- trim_whitespace: removes leading whitespace characters from a string in-place.
- str2array: Parses a string containing comma-separated values enclosed in brackets and converts them to a double array.
*/

void trim_whitespace(char* s) {
  const char* d = s;
  do {
    while (*d == ' ')
      ++d;
  } while ((*s++ = *d++));
  *s = '\0';
}


void str2array(char *tmps2, double *array){

  // overwrite newline character with string terminator
  char *newline = strchr( tmps2, '\n' );
  if ( newline )
    *newline = 0;

  char* p;
  int n = 0;

  p = strtok(tmps2,"[,]");
  while (p != NULL){
    array[n] = atof(p);
    fprintf (stderr, "array[%d] = %g \n", n, array[n]);
    p = strtok(NULL, ",");
    n += 1;
  }
}

/**
   Params global variable
*/

Array * params;

typedef struct {
  char * name;
  void * ptr;
  char * type;
  int len;
} ParamItem;


void add_param(char * name, void * ptr,  char * type, int len = 0)
{
  ParamItem p;

  p.name = strdup(name);
  p.ptr = ptr;
  p.type = strdup(type);
  p.len = len;

  array_append (params, &p, sizeof (ParamItem));

}

// I prefer not to put these 2 in an event because the declaration must be done
// before run(). (same rule for cleanup for simplicity)
void init_param()
{
params = array_new();
}

void cleanup_param()
{
  ParamItem * d2 = params->p;
  for (int i = 0; i < params->len/sizeof(ParamItem); i++, d2++) { // loop over parameters
    free(d2->name);
    free(d2->type);
  }
  array_free (params);
}


void read_params(char* path2file)
{
  FILE * fp;
  if ((fp = fopen(path2file, "rt"))) {
    char tempbuff[300];
    while(fgets(tempbuff,300,fp)) { // loop over file lines
      trim_whitespace(tempbuff);
      char* tmps1 = strtok(tempbuff, "=");
      char* tmps2 = strtok(NULL, "=");

  ParamItem * d2 = params->p;
  for (int i = 0; i < params->len/sizeof(ParamItem); i++, d2++) { // loop over parameters
    if (strcmp(d2->name, tmps1) == 0) {
      //check for type and assign value
    if (strcmp(d2->type, "int") == 0) {
      *( (int*) d2->ptr) = atoi(tmps2);
      fprintf (stderr, "scan param %s: %d\n", d2->name, *( (int*) d2->ptr));}
    else if (strcmp(d2->type, "double") == 0){
      *( (double*) d2->ptr) = atof(tmps2);
      fprintf (stderr, "scan param %s: %e\n", d2->name, *( (double*) d2->ptr));}
    else if (strcmp(d2->type, "array") == 0){
      fprintf (stderr, "scan param %s:\n", d2->name);
      str2array(tmps2, (double*) d2->ptr);
      }
    }
    }

    }
    fclose(fp);
  } else {
    fprintf(stdout, "file %s not found\n", path2file);
    exit(0);
  }
}

/**
   Create output directory and copy input parameter file for backup
*/
void create_outdir()
{
  if (pid() == 0) {
    for (int i=1; i<10000; i++) {
      sprintf(dpath, "outdir_%04d/", i);
      if (mkdir(dpath, 0777) == 0) {
        fprintf(stdout,"Writing output in %s\n",dpath);
        break;
      }
    }
  }
@if _MPI
  MPI_Bcast(&dpath, 80, MPI_CHAR, 0, MPI_COMM_WORLD);
@endif
}

void backup_config(char* path2file)
{
  fprintf(stdout, "Backup config\n");
  if (pid() == 0) {
    char ch;
    char name[90];
    sprintf (name,"%sparams.in", dpath);
    FILE * source = fopen(path2file, "r");
    FILE * target = fopen(name, "w");
    while ((ch = fgetc(source)) != EOF)
      fputc(ch, target);
    fclose(source);
    fclose(target);
  }
}

/**
 Copy file
 from: https://stackoverflow.com/questions/29079011/copy-file-function-in-c
 */
void backup_file(char FileSource[])
{
  if (pid() == 0) {
    char FileDestination[100];
    sprintf (FileDestination,"%s%s", dpath, FileSource);

    char    c[4096]; // or any other constant you like
    FILE    *stream_R = fopen(FileSource, "r");
    FILE    *stream_W = fopen(FileDestination, "w");   //create and write to file

    while (!feof(stream_R)) {
      size_t bytes = fread(c, 1, sizeof(c), stream_R);
      if (bytes) {
        fwrite(c, 1, bytes, stream_W);
      }
    }

    //close streams
    fclose(stream_R);
    fclose(stream_W);
  }
}
/**
## Test
* [Skeleton file for params handling](../tests/test_params.c)

*/
