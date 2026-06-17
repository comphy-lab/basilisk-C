#include <stdint.h>

typedef struct {
  uint32_t s;
} Adler32Hash;

static
inline void a32_hash_init (Adler32Hash * hash)
{
  hash->s = 0;
}

static
inline void a32_hash_add (Adler32Hash * hash, const void * data, size_t size)
{
  const uint8_t * buffer = (const uint8_t*) data;
  for (size_t n = 0; n < size; n++, buffer++)
    hash->s = *buffer + (hash->s << 6) + (hash->s << 16) - hash->s;
}

static
inline uint32_t a32_hash (const Adler32Hash * hash)
{
  return hash->s;
}
