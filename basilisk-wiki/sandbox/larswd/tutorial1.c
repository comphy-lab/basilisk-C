/**
# A small crash course in Basilisk C (multilayer version)
Basilisk C is a programming language ( a dialect of C if you will) made for solving partial differential equations. The software is useful both for coding your own custom-made solvers as well as importing pre-existing solvers. This tutorial will focus on one of the many preprogrammed PDE solvers in Basilisk C, the multilayer solver. But before we begin on the multilayer solver, we will start by first giving a quick tutorial of C programming before moving on to a quick tutorial of Basilisk C. 

## C programming in a shellnut
This section is written for users who are familiar with programming languages such as python. Hence, if you are familiar with C programming you can skip ahead to the [next section](tutorial3.c).

The C programming language is a quite old and well known (perhaps even notorious) programming language. Unlike more "modern" programming languages such as python, matlab, or R, C is a compiled language with static typing. A compiled language means that you do not run the code file you have written directly. Instead, you must first use a compiler, a translation program, to convert the script to an executable file which you can then run. In linux operating systems, (and the windows subsystem for linux) you have access to the ```gcc``` compiler. To see that this is in order and working, write 

```bash
gcc
```
in your terminal. You should see something like:
```bash
$ gcc
gcc: fatal error: no input files
compilation terminated.
```
Which shows that gcc is installed on your system. 
*/