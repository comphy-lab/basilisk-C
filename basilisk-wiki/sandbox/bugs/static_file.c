/**
# Declaring and assignment of static FILE crashes qcc

see also [here on the forum](https://groups.google.com/g/basilisk-fr/c/kfr0SwsPtHg)
*/
int main() {
  static FILE * fp = fopen ("file", "w");
}