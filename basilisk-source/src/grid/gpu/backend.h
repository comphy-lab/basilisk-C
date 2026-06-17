// fixme: for the moment only 'const int' are considered, this could be generalised
#define IS_EXTERNAL_CONSTANT(g) ((g)->constant && (g)->type == sym_INT && !(g)->data)
#define EXTERNAL_NAME(g) (g)->global == 2 ? "_loc_" : "", (g)->name, (g)->reduct ? "_in_" : ""

static char * str_append_array (char * dst, const char * list[])
{
  int empty = (dst == NULL);
  int len = empty ? 0 : strlen (dst);
  for (const char ** s = list; *s != NULL; s++)
    len += strlen (*s);
  dst = (char *) sysrealloc (dst, len + 1);
  if (empty) dst[0] = '\0';
  for (const char ** s = list; *s != NULL; s++)
    strcat (dst, *s);
  return dst;
}

#define str_append(dst, ...) str_append_array (dst, (const char *[]){__VA_ARGS__, NULL})
#define xstr(a) str(a)
#define str(a) #a

#include <khash.h>

typedef struct _Shader Shader;
typedef struct _GPUData GPUData;

KHASH_MAP_INIT_INT(INT, Shader *)

void gpu_free_solver (void);
void free_shader (Shader * s);
Shader * load_normal_shader (const char * fs, const char * func, const char * file, int line);
void gpu_free_context (GPUData * data);
bool gpu_init_context (GPUData ** data);
void realloc_ssbo (size_t field_size);
void gpu_synchronize();

typedef enum {
  GPU_READ, GPU_WRITE
} SyncMode;

void gpu_cpu_sync_scalar (int i, int block, char * data, size_t field_size, SyncMode mode);
void reset_scalar (int i, int block, size_t field_size, double val);

void finalize_shader (Shader * s, External * externals, External * merged,
                      unsigned ng[2], unsigned nwg[2]);
void post_setup_shader (Shader * shader, External * externals);
int run_shader (const Shader * shader, const RegionParameters * region);
double gpu_reduction (size_t offset, const char op, const RegionParameters * region,
                      GPUData * data, size_t nb);
char * gpu_errors (const char * errors, const char * source, char * fout,
                   const char * lang);

void printWorkGroupsCapabilities();
void gpu_limits (FILE * fp);

typedef struct {
  bool fragment_shader;
  int current_shader, nssbo;
  size_t max_ssbo_size, current_size;  
} GPUContext_t;

extern GPUContext_t GPUContext;

#undef device_synchronize
#define device_synchronize() gpu_synchronize()
