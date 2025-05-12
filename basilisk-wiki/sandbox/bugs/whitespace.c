/**
# Whitespace in the syntax

qcc forgets to add whitespaces if the user doesn't when declaring boundary conditions.
*/

scalar f[];
f[left]=0;

int main() { 
  ;
}