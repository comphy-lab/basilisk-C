/**
# Radix sort

Modified from Geeks for Geeks:  
[https://www.geeksforgeeks.org/radix-sort/]()

Sorting positive integers, keeping track of their original index.

A linear helper algorithm:  
*/
void countSort (int n, int arr[n], int ind[n], int exp) { 
  int output[n], out_ind[n]; 
    int i, count[10] = {0}; 
    for (i = 0; i < n; i++) 
        count[(arr[i] / exp) % 10]++; 
    for (i = 1; i < 10; i++) 
        count[i] += count[i - 1]; 
    for (i = n - 1; i >= 0; i--) {
      out_ind[count[(arr[i] / exp) % 10] - 1] = ind[i];
      output[count[(arr[i] / exp) % 10] - 1]  = arr[i]; 
      count[(arr[i] / exp) % 10]--; 
    } 
    for (i = 0; i < n; i++) { 
      ind[i] = out_ind[i];
      arr[i] = output[i];
    }
}
/**
## User interface

Input:  

n   , number of elements in array  
arr , the array to be sorted (input and output)  
Ind , index translation array (output)  
*/
void radixsort (int n, int arr[n], int ind[n]) {
  int m = arr[0]; 
  for (int i = 0; i < n; i++) {
    if (arr[i] > m) 
      m = arr[i];
    ind[i] = i;
  }
  for (int exp = 1; m/exp > 0; exp *= 10)
    countSort(n, arr, ind, exp);
}
