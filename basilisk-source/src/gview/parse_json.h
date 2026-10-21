#include <stdio.h>
#include <string.h>

int parse_entry (char *s, char **name, char **keys, char **values, int max_pairs)
{
  int n = 0;

  /* Name */
  if (*s++ != '"')
    return -1;
  *name = s;
  s = strchr (s, '"');
  if (!s)
    return -1;
  *s++ = '\0';

  /* Find first key */
  s = strchr (s, '{');
  if (!s)
    return -1;
  s++;

  while (n < max_pairs) {
    /* Key */
    s = strchr (s, '"');
    if (!s || s[1] == '}')
      break;

    keys[n] = ++s;
    s = strchr (s, '"');
    if (!s)
      return -1;
    *s++ = '\0';

    /* Value */
    s = strchr (s, ':');
    if (!s)
      return -1;
    s++;

    while (*s == ' ' || *s == '\t')
      s++;

    if (*s == '"') {
      values[n] = ++s;
      s = strchr (s, '"');
      if (!s)
        return -1;
      *s++ = '\0';
    }
    else {
      values[n] = s;
      s = strpbrk (s, ",}");
      if (!s)
        return -1;
      *s++ = '\0';
    }

    /* Remove trailing whitespace from unquoted value */
    char *end = values[n] + strlen (values[n]);
    while (end > values[n] &&
           (end[-1] == ' ' || end[-1] == '\t'))
      *--end = '\0';

    n++;
  }

  return n;
}
