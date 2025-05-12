int main() {
  uint32_t abc = 123;
  printf("a = %" PRIu32 "\n", abc);
  //printf("a = %u\n", abc);
}

/*
$ qcc main.c -lm
main.c:3: warning: Basilisk C parse error near `printf("a = %" PRIu32 "\n", abc)'
main.c: In function ‘main’:
main.c:3:17: error: expected ‘)’ before ‘PRIu32’
main.c:3:10: warning: spurious trailing ‘%’ in format [-Wformat=]
*/
