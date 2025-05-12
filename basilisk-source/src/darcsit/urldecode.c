#include <stdio.h>
#include <ctype.h>

int main()
{
  int c = getchar(), c1 = getchar(), c2 = getchar();
  while (c != EOF) {
    // fprintf (stderr, "%c|%c|%c\n", c, c1, c2);
    if (c == '%' && isxdigit(c1) && isxdigit(c2)) {
      int a = c1, b = c2;
      if (a >= 'a')
	a -= 'a'-'A';
      if (a >= 'A')
	a -= ('A' - 10);
      else
	a -= '0';
      if (b >= 'a')
	b -= 'a'-'A';
      if (b >= 'A')
	b -= ('A' - 10);
      else
	b -= '0';
      putchar (16*a + b);
      c = getchar(), c1 = getchar(), c2 = getchar();
    } else if (c == '+') {
      putchar (' ');
      c = c1, c1 = c2, c2 = getchar();
    } else {
      putchar (c);
      c = c1, c1 = c2, c2 = getchar();
    }
  }
  return 0;
}
