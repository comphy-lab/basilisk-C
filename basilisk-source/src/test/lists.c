// Check that lists work properly
// some compilers are picky about how lists of structures can be
// allocated

scalar sa[];
vector sb[];

scalar * list2 = {sa};
vector * list3 = {sb};
scalar * list4 = {sa,sb.x,sb.y};

tensor tc[];
scalar * list5 = {sa,tc.x.x,tc.x.y,tc.y.x,tc.y.y};
vector * list6 = {sb,tc.x,tc.y};
scalar * list7 = {sa,tc.x,tc.y};
scalar * list8 = {sa,tc};
scalar * list9 = {sa,sb};

int main()
{
  for (scalar s in list4)
    fprintf (stderr, "%d\n", s.i);
  for (scalar s in list5)
    fprintf (stderr, "%d\n", s.i);
  for (vector v in list6)
    fprintf (stderr, "%d %d\n", v.x.i, v.y.i);
  for (scalar s in list7)
    fprintf (stderr, "%d\n", s.i);
  for (scalar s in list8)
    fprintf (stderr, "%d\n", s.i);
  for (scalar s in list9)
    fprintf (stderr, "%d\n", s.i);
  
  scalar a = {1}, b = {2}, c = {3};
  vector d = {{4}, {5}};
  scalar * list = {a,b,c,d};
  for (scalar s in list)
    fprintf (stderr, "%d\n", s.i);

  tensor e = {{{6},{7}},{{8},{9}}};
  vector * list1 = {d,e};
  for (vector v in list1)
    fprintf (stderr, "%d %d\n", v.x.i, v.y.i);

  tensor * list2 = {e};
  for (tensor t in list2)
    fprintf (stderr, "%d %d %d %d\n", t.x.x.i, t.x.y.i, t.y.x.i, t.y.y.i);

  scalar * list3 = {a,d.x,d.y};
  for (scalar s in list3)
    fprintf (stderr, "%d\n", s.i);
}
