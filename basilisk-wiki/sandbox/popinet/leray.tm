<TeXmacs|1.99.2>

<style|generic>

<\body>
  Start from the incompressible Euler equations

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<partial\><rsub|t>\<b-u\>+\<b-u\>\<cdot\><math-bf|\<nabla\>>\<b-u\>>|<cell|=>|<cell|-<math-bf|\<nabla\>>p>>|<row|<cell|<math-bf|\<nabla\>>\<cdot\>\<b-u\>>|<cell|=>|<cell|0>>>>
  </eqnarray*>

  Assume a self-similar solution

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<b-u\><around*|(|\<b-r\>,t|)>>|<cell|=>|<cell|<around*|(|t<rsup|\<star\>>-t|)><rsup|-\<alpha\>>*\<b-U\><around*|(|\<b-r\>*<around*|(|t<rsup|\<star\>>-t|)><rsup|-\<beta\>>|)>>>>>
  </eqnarray*>

  For <math|\<alpha\>=\<beta\>=1/2>, this gives equations (3) and (4) of
  Pomeau (i.e. Leray without viscosity)

  <\eqnarray*>
    <tformat|<table|<row|<cell|<frac|1|2>*<around*|(|\<b-U\>+\<b-R\>\<cdot\><math-bf|\<nabla\>>\<b-U\>|)>+\<b-U\>\<cdot\><math-bf|\<nabla\>>\<b-U\>>|<cell|=>|<cell|-<math-bf|\<nabla\>>P>>|<row|<cell|<math-bf|\<nabla\>>\<cdot\>\<b-U\>>|<cell|=>|<cell|0>>>>
  </eqnarray*>

  This can be interpreted as a stationary solution of the equations

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<partial\><rsub|t<rsup|\<star\>>>\<b-U\>+<around*|(|\<b-U\>+<frac|1|2>*\<b-R\>|)>\<cdot\><math-bf|\<nabla\>>\<b-U\>>|<cell|=>|<cell|-<math-bf|\<nabla\>>P-<frac|1|2>*\<b-U\>>>|<row|<cell|<math-bf|\<nabla\>>\<cdot\>\<b-U\>>|<cell|=>|<cell|0>>>>
  </eqnarray*>

  with <math|<rsup|>t<rsup|\<star\>>> a pseudo-time, i.e. the incompressible
  Euler equations with linear damping (<math|><math|-<frac|1|2>*\<b-U\>>) and
  an added radial expansion velocity <math|<frac|1|2>*\<b-R\>>. This
  evolution equation can be discretised numerically using standard schemes.
  Unstable stationary solutions could in principle be obtained through
  stabilisation by relaxation.

  A non-trivial stationary solution exists only when the radial expansion
  term is balanced by some (compressive) term. This compression cannot come
  from the velocity field itself (since it is divergence-free). I cannot see
  how it could come from the linear damping term either. So it seems that
  non-trivial stationary solutions do not exist. Is this correct?

  If one considers a stationary solution <math|\<b-U\><rsub|0>> of the
  incompressible Euler equations, we have

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<b-U\><rsub|0>\<cdot\><math-bf|\<nabla\>>\<b-U\><rsub|0>>|<cell|=>|<cell|-<math-bf|\<nabla\>>P>>|<row|<cell|<math-bf|\<nabla\>>\<cdot\>\<b-U\><rsub|0>>|<cell|=>|<cell|0>>>>
  </eqnarray*>

  Using this as initial condition gives at <math|t=0>

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<partial\><rsub|t>\<b-U\><rsub|0>+<frac|1|2>*\<b-R\>\<cdot\><math-bf|\<nabla\>>\<b-U\><rsub|0>>|<cell|=>|<cell|-<frac|1|2>*\<b-U\><rsub|0>>>|<row|<cell|<math-bf|\<nabla\>>\<cdot\>\<b-U\><rsub|0>>|<cell|=>|<cell|0>>>>
  </eqnarray*>

  If <math|\<b-U\><around*|(|\<b-R\>|)>> is a solution then
  <math|\<lambda\>*\<b-U\><around*|(|\<lambda\>*\<b-R\>|)>> is also a
  solution.

  The pressure is given by

  <\eqnarray*>
    <tformat|<table|<row|<cell|<frac|1|2>*<math-bf|\<nabla\>>\<cdot\><around*|(|\<b-R\>\<cdot\><math-bf|\<nabla\>>\<b-U\>|)>+<math-bf|\<nabla\>>\<cdot\><around*|(|\<b-U\>\<cdot\><math-bf|\<nabla\>>\<b-U\>|)>>|<cell|=>|<cell|-<math-bf|\<nabla\>><rsup|2>P>>|<row|<cell|<math-bf|\<nabla\>>\<cdot\><around*|(|\<b-U\>\<cdot\><math-bf|\<nabla\>>\<b-U\>|)>>|<cell|=>|<cell|-<math-bf|\<nabla\>><rsup|2>P>>>>
  </eqnarray*>

  Solution procedure 1

  <\enumerate-numeric>
    <item>Compute <math|\<b-U\><rsup|\<star\>>>

    <\equation*>
      \<b-U\><rsup|\<star\>>=-<around*|(|2*\<b-U\><rsup|n>+\<b-R\>|)>\<cdot\><math-bf|\<nabla\>>\<b-U\><rsup|n>
    </equation*>

    <item>Projection: <math|\<b-U\><rprime|'>=\<b-U\><rsup|\<star\>>-<math-bf|\<nabla\>>P<rprime|'>>
    with <math|><math|<math-bf|\<nabla\>>\<cdot\>\<b-U\><rprime|'>=0> and

    <\equation*>
      <math-bf|\<nabla\>><rsup|2>P<rprime|'>=<math-bf|\<nabla\>>\<cdot\>\<b-U\><rsup|\<star\>>
    </equation*>

    <item>Dilation map

    <\equation*>
      \<b-U\><rsup|n+1>=\<lambda\>*\<b-U\><rprime|'><around*|(|\<lambda\>*\<b-R\>|)>
    </equation*>
  </enumerate-numeric>

  Solution procedure 2 (from Pomeau)

  <\enumerate-numeric>
    <item>Compute <math|\<b-K\>>

    <\equation*>
      \<b-K\><rsup|n>=\<b-U\><rsup|n>\<cdot\><math-bf|\<nabla\>>\<b-U\><rsup|n>
    </equation*>

    <item>Compute <math|P>

    <\equation*>
      <math-bf|\<nabla\>><rsup|2>P<rsup|n>=-<math-bf|\<nabla\>>\<cdot\>\<b-K\><rsup|n>
    </equation*>

    <item>Compute <math|\<b-V\>>

    <\equation*>
      \<b-V\><rsup|n>=-2*<around*|(|<math-bf|\<nabla\>>P<rsup|n>+\<b-K\><rsup|n>|)>
    </equation*>

    <item>Compute <math|<wide|\<b-U\>|~>>

    <\equation*>
      <wide|\<b-U\>|~>+\<b-R\>\<cdot\><math-bf|\<nabla\>><wide|\<b-U\>|~>=\<b-V\><rsup|n>
    </equation*>

    which gives

    <\equation*>
      <wide|\<b-U\>|~><around*|(|\<b-R\>|)>=<frac|1|R>*<big|int><rsup|R><rsub|0>\<b-V\><rsup|n><around*|(|R<rprime|'>,<wide|\<b-R\>|^>|)>*d
      R<rprime|'>
    </equation*>

    <item>Dilation map

    <\equation*>
      \<b-U\><rsup|n+1>=\<lambda\>*<wide|\<b-U\>|~><around*|(|\<lambda\>*\<b-R\>|)>
    </equation*>
  </enumerate-numeric>

  Linear first-order ODE

  <\eqnarray*>
    <tformat|<table|<row|<cell|u+x*u<rprime|'>>|<cell|=>|<cell|a>>|<row|<cell|u+<around*|(|x*u|)><rprime|'>-u>|<cell|=>|<cell|a>>|<row|<cell|<around*|(|x*u|)><rprime|'>>|<cell|=>|<cell|a>>|<row|<cell|u>|<cell|=>|<cell|<frac|1|x>*<big|int><rsup|x><rsub|0>a<around*|(|x<rprime|'>|)>*d
    x<rprime|'>>>>>
  </eqnarray*>

  Example:

  <\eqnarray*>
    <tformat|<table|<row|<cell|a>|<cell|=>|<cell|cos<around*|(|x|)>>>|<row|<cell|u>|<cell|=>|<cell|<frac|sin<around*|(|x|)>|x>>>>>
  </eqnarray*>

  \;

  \;

  \;
</body>

<initial|<\collection>
</collection>>