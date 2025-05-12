/**
# Small typos crash qcc

A usefull compiler error message is expected but qcc crashes instead. Making debugging  harder than is should be.
*/

int main () {
  double a = 5 // missing `;`
  }

/**
A seemingly similar issue with this code:

~~~c
int main () {
  { 
    double a = 5;
  }
// Missing `}`
~~~
*/