<TeXmacs|1.99.2>

<style|generic>

<\body>
  <section|Steady-state force balance>

  Consider a domain <math|\<Omega\>> containing an obstacle <math|\<Gamma\>>,
  at steady state we have the balance

  <\equation*>
    <big|int><rsub|\<Omega\>>-<math-bf|\<nabla\>>p+<math-bf|\<nabla\>>\<cdot\><around*|(|2*\<mu\>*\<b-D\>|)>=0
  </equation*>

  which can also be written

  <\equation*>
    \<b-F\><rsub|\<Gamma\>>=<big|int><rsub|\<partial\>\<Gamma\>><around*|(||\<nobracket\>>-p*\<b-I\>+<around*|\<nobracket\>|2*\<mu\>*\<b-D\>|)>\<cdot\>\<b-n\>*d\<partial\>\<Gamma\>=-<big|int><rsub|\<partial\>\<Omega\>><around*|(||\<nobracket\>>-p*\<b-I\>+<around*|\<nobracket\>|2*\<mu\>*\<b-D\>|)>\<cdot\>\<b-n\>*d\<partial\>\<Omega\>
  </equation*>

  i.e. the force on the obstacle can be obtained by integrating the stresses
  on the (much simpler) boundary of <math|\<Omega\>>.

  From Choi et al. (2007)

  <\eqnarray*>
    <tformat|<table|<row|<cell|F<rsub|i>>|<cell|=>|<cell|<big|int><rsub|V>\<partial\><rsub|t><around*|(|\<rho\>*u<rsub|i>|)>*d
    V+<big|int><rsub|S<rsub|o>><around*|(|\<rho\>*u<rsub|i>*u<rsub|j>+p*\<delta\><rsub|i
    j>-\<tau\><rsub|i j>|)>*n<rsub|j>*d s<rsub|o>>>>>
  </eqnarray*>

  This approach was used by:

  J.-I. Choi et al. / Journal of Computational Physics 224 (2007) 757\U784

  E. Balaras, Modeling complex boundaries using an external force field on
  fixed Cartesian grids in large-eddy simulations, Computers and Fluids 33
  (2004) 375\U404.

  <section|Direct integration>

  <\equation*>
    \<b-D\>=<frac|1|2>*<around*|(|<math-bf|\<nabla\>>\<b-u\>+<math-bf|\<nabla\>><rsup|T>\<b-u\>|)>=<frac|1|2>*<matrix|<tformat|<table|<row|<cell|2*\<partial\><rsub|x>u>|<cell|\<partial\><rsub|x>v+\<partial\><rsub|y>u>>|<row|<cell|\<partial\><rsub|x>v+\<partial\><rsub|y>u>|<cell|2*\<partial\><rsub|y>v>>>>>
  </equation*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<b-F\><rsub|\<Gamma\>>>|<cell|=>|<cell|-<big|int><rsub|\<Gamma\>><matrix|<tformat|<table|<row|<cell|<around*|(|-p+2*\<mu\>*\<partial\><rsub|x>u|)>*n<rsub|x>+\<mu\>*<around*|(|\<partial\><rsub|x>v+\<partial\><rsub|y>u|)>*n<rsub|y>>>|<row|<cell|<around*|(|-p+2*\<mu\>*\<partial\><rsub|y>v|)>*n<rsub|y>+\<mu\>*<around*|(|\<partial\><rsub|x>v+\<partial\><rsub|y>u|)>*n<rsub|x>>>>>>>>|<row|<cell|>|<cell|=>|<cell|-<big|int><rsub|\<Gamma\>><matrix|<tformat|<table|<row|<cell|-p*n<rsub|x>+2*\<mu\>*<math-bf|\<nabla\>>u\<cdot\>\<b-n\>+\<mu\>*<around*|(|\<partial\><rsub|x>v-\<partial\><rsub|y>u|)>*n<rsub|y>>>|<row|<cell|-p*n<rsub|y>+2*\<mu\>*<math-bf|\<nabla\>>v\<cdot\>\<b-n\>-\<mu\>*<around*|(|\<partial\><rsub|x>v-\<partial\><rsub|y>u|)>*n<rsub|x>>>>>>>>>>
  </eqnarray*>

  If we assume that <math|\<b-u\>> is constant on the boundary, then

  <\equation*>
    <math-bf|\<nabla\>>\<b-u\>\<cdot\>\<b-t\>=<math-bf|0>
  </equation*>

  with <math|\<b-t\>> the unit tangent vector to the boundary. We thus have
  the relations

  <\eqnarray*>
    <tformat|<table|<row|<cell|<math-bf|\<nabla\>>\<b-u\>>|<cell|=>|<cell|<around*|(|<math-bf|\<nabla\>>\<b-u\>\<cdot\>\<b-n\>|)>*\<b-n\>+<around*|(|<math-bf|\<nabla\>>\<b-u\>\<cdot\>\<b-t\>|)>*\<b-t\>=<around*|(|<math-bf|\<nabla\>>\<b-u\>\<cdot\>\<b-n\>|)>*\<b-n\>>>|<row|<cell|\<b-D\>=<frac|1|2>*<around*|(|<math-bf|\<nabla\>>\<b-u\>+<math-bf|\<nabla\>><rsup|T>\<b-u\>|)>>|<cell|=>|<cell|<frac|1|2>*<matrix|<tformat|<table|<row|<cell|2*<around*|(|<math-bf|\<nabla\>>u\<cdot\>\<b-n\>|)>*n<rsub|x>>|<cell|<around*|(|<math-bf|\<nabla\>>u\<cdot\>\<b-n\>|)>*n<rsub|y>+<around*|(|<math-bf|\<nabla\>>v\<cdot\>\<b-n\>|)>*n<rsub|x>>>|<row|<cell|<around*|(|<math-bf|\<nabla\>>u\<cdot\>\<b-n\>|)>*n<rsub|y>+<around*|(|<math-bf|\<nabla\>>v\<cdot\>\<b-n\>|)>*n<rsub|x>>|<cell|2*<around*|(|<math-bf|\<nabla\>>v\<cdot\>\<b-n\>|)>*n<rsub|y>>>>>>>>|<row|<cell|\<b-F\><rsub|\<Gamma\>>>|<cell|=>|<cell|-<big|int><rsub|\<Gamma\>><matrix|<tformat|<table|<row|<cell|<around*|[|-p+2*\<mu\>*<around*|(|<math-bf|\<nabla\>>u\<cdot\>\<b-n\>|)>*n<rsub|x>|]>*n<rsub|x>+\<mu\>*<around*|[|<around*|(|<math-bf|\<nabla\>>u\<cdot\>\<b-n\>|)>*n<rsub|y>+<around*|(|<math-bf|\<nabla\>>v\<cdot\>\<b-n\>|)>*n<rsub|x>|]>*n<rsub|y>>>|<row|<cell|<around*|[|-p+2*\<mu\>*<around*|(|<math-bf|\<nabla\>>v\<cdot\>\<b-n\>|)>*n<rsub|y>|]>*n<rsub|y>+\<mu\>*<around*|[|<around*|(|<math-bf|\<nabla\>>u\<cdot\>\<b-n\>|)>*n<rsub|y>+<around*|(|<math-bf|\<nabla\>>v\<cdot\>\<b-n\>|)>*n<rsub|x>|]>*n<rsub|x>>>>>>>>|<row|<cell|>|<cell|=>|<cell|-<big|int><rsub|\<Gamma\>><matrix|<tformat|<table|<row|<cell|-p*n<rsub|x>+\<mu\>*<around*|[|<around*|(|<math-bf|\<nabla\>>u\<cdot\>\<b-n\>|)>*<around*|(|n<rsup|2><rsub|x>+1|)>+<around*|(|<math-bf|\<nabla\>>v\<cdot\>\<b-n\>|)>*n<rsub|x>*n<rsub|y>|]>>>|<row|<cell|-p*n<rsub|y>+\<mu\>*<around*|[|<around*|(|<math-bf|\<nabla\>>v\<cdot\>\<b-n\>|)>*<around*|(|n<rsup|2><rsub|y>+1|)>+<around*|(|<math-bf|\<nabla\>>u\<cdot\>\<b-n\>|)>*n<rsub|x>*n<rsub|y>|]>>>>>>>>>>
  </eqnarray*>

  The vorticity is

  <\equation*>
    \<omega\>=\<partial\><rsub|x>v-\<partial\><rsub|y>u=<around*|(|<math-bf|\<nabla\>>v\<cdot\>\<b-n\>|)>*n<rsub|x>-<around*|(|<math-bf|\<nabla\>>u\<cdot\>\<b-n\>|)>*n<rsub|y>
  </equation*>

  \;
</body>

<initial|<\collection>
</collection>>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|?>>
    <associate|auto-2|<tuple|2|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Steady-state
      force balance> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Direct
      integration> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>