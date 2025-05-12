<TeXmacs|1.0.7.19>

<style|article>

<\body>
  <section|Balanced Saint-Venant scheme>

  <subsection|Uniform case>

  From Audusse et al, 2005, the one-dimensional scheme in Cartesian
  coordinates can be written

  <\equation>
    \<Delta\>x<rsub|i>*d<rsub|t>U<rsub|i>+\<cal-F\><rsub|l>-\<cal-F\><rsub|r>=0,<label|audusse>
  </equation>

  with left and right numerical fluxes

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<cal-F\><rsub|l>>|<cell|=>|<cell|\<cal-F\><around|(|U<rsub|i+1/2->,U<rsub|i+1/2+>|)>+<frac|g|2>*<matrix|<tformat|<table|<row|<cell|0>>|<row|<cell|h<rsup|2><rsub|i,r>-h<rsup|2><rsub|i+1/2->+<around*|(|h<rsub|i,r>+h<rsub|i>|)>*<around*|(|z<rsub|i,r>-z<rsub|i>|)>>>>>>,>>|<row|<cell|\<cal-F\><rsub|r>>|<cell|=>|<cell|\<cal-F\><around|(|U<rsub|i-1/2->,U<rsub|i-1/2+>|)>+<frac|g|2>*<matrix|<tformat|<table|<row|<cell|0>>|<row|<cell|h<rsup|2><rsub|i,l>-h<rsup|2><rsub|i-1/2+>+<around*|(|h<rsub|i,l>+h<rsub|i>|)>*<around*|(|z<rsub|i,l>-z<rsub|i>|)>>>>>>,>>>>
  </eqnarray*>

  The lake-at-rest steady state is defined by <math|u<rsub|i>=0> and
  <math|h<rsub|i,l>+z<rsub|i,l>=h<rsub|i,r>+z<rsub|i,r>=h<rsub|i>+z<rsub|i>=H>
  for all <math|i>. In this case

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<cal-F\><rsub|l>>|<cell|=>|<cell|<frac|g|2>*<matrix|<tformat|<table|<row|<cell|0>>|<row|<cell|*\h<rsup|2><rsub|i+1/2->>>>>>+<frac|g|2>*<matrix|<tformat|<table|<row|<cell|0>>|<row|<cell|\h<rsup|2><rsub|i,r>-\h<rsup|2><rsub|i+1/2->+<around*|(|\h<rsub|i,r>+\h<rsub|i>|)>*<around*|(|\h<rsub|i>-\h<rsub|i,r>|)>>>>>>,>>|<row|<cell|\<cal-F\><rsub|r>>|<cell|=>|<cell|<frac|g|2>*<matrix|<tformat|<table|<row|<cell|0>>|<row|<cell|\h<rsup|2><rsub|i-1/2+>>>>>>+<frac|g|2>*<matrix|<tformat|<table|<row|<cell|0>>|<row|<cell|*\h<rsup|2><rsub|i,l>-\h<rsup|2><rsub|i-1/2+>+<around*|(|\h<rsub|i,l>+\h<rsub|i>|)>*<around*|(|\h<rsub|i>-\h<rsub|i,l>|)>>>>>>,>>>>
  </eqnarray*>

  so that

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<cal-F\><rsub|l>-\<cal-F\><rsub|r>>|<cell|=>|<cell|<frac|g|2>*<matrix|<tformat|<table|<row|<cell|0>>|<row|<cell|\h<rsup|2><rsub|i,r>+\h<rsup|2><rsub|i>-\h<rsup|2><rsub|i,r>-*\h<rsup|2><rsub|i,l>+\h<rsup|2><rsub|i,l>-\h<rsup|2><rsub|i>>>>>>=<matrix|<tformat|<table|<row|<cell|0>>|<row|<cell|*0>>>>>>>>>
  </eqnarray*>

  <subsection|Fine/coarse case>

  From Audusse et al, 2005, the one-dimensional scheme in Cartesian
  coordinates can be written

  <\equation>
    \<Delta\>x<rsub|i>*d<rsub|t>U<rsub|i>+\<cal-F\><rsub|l>-\<cal-F\><rsub|r>=0,<label|audusse>
  </equation>

  with left and right numerical fluxes

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<cal-F\><rsub|l>>|<cell|=>|<cell|\<cal-F\><around|(|U<rsub|i+1/2->,U<rsub|i+1/2+>|)>+<frac|g|2>*<matrix|<tformat|<table|<row|<cell|0>>|<row|<cell|\h<rsup|2><rsub|i,r>-\h<rsup|2><rsub|i+1/2->+<around*|(|\h<rsub|i,r>+\h<rsub|i,f>|)>*<around*|(|z<rsub|i,r>-z<rsub|i,f>|)>>>>>>,>>|<row|<cell|\<cal-F\><rsub|r>>|<cell|=>|<cell|\<cal-F\><around|(|U<rsub|i-1/2->,U<rsub|i-1/2+>|)>+<frac|g|2>*<matrix|<tformat|<table|<row|<cell|0>>|<row|<cell|\h<rsup|2><rsub|i,l>-\h<rsup|2><rsub|i-1/2+>+<around*|(|\h<rsub|i,l>+\h<rsub|i,c>|)>*<around*|(|z<rsub|i,l>-z<rsub|i,c>|)>>>>>>,>>>>
  </eqnarray*>

  The lake-at-rest steady state is defined by <math|u<rsub|i>=0> and
  <math|\h<rsub|i,l>+z<rsub|i,l>=\h<rsub|i,r>+z<rsub|i,r>=\h<rsub|i,f>+z<rsub|i,f>=\h<rsub|i,c>+z<rsub|i,c>=H>
  for all <math|i>. In this case

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<cal-F\><rsub|l>>|<cell|=>|<cell|<frac|g|2>*<matrix|<tformat|<table|<row|<cell|0>>|<row|<cell|*\h<rsup|2><rsub|i+1/2->>>>>>+<frac|g|2>*<matrix|<tformat|<table|<row|<cell|0>>|<row|<cell|\h<rsup|2><rsub|i,r>-\h<rsup|2><rsub|i+1/2->+<around*|(|\h<rsub|i,r>+\h<rsub|i,f>|)>*<around*|(|\h<rsub|i,f>-\h<rsub|i,r>|)>>>>>>,>>|<row|<cell|\<cal-F\><rsub|r>>|<cell|=>|<cell|<frac|g|2>*<matrix|<tformat|<table|<row|<cell|0>>|<row|<cell|\h<rsup|2><rsub|i-1/2+>>>>>>+<frac|g|2>*<matrix|<tformat|<table|<row|<cell|0>>|<row|<cell|*\h<rsup|2><rsub|i,l>-\h<rsup|2><rsub|i-1/2+>+<around*|(|\h<rsub|i,l>+\h<rsub|i,c>|)>*<around*|(|\h<rsub|i,c>-\h<rsub|i,l>|)>>>>>>,>>>>
  </eqnarray*>

  so that

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<cal-F\><rsub|l>-\<cal-F\><rsub|r>>|<cell|=>|<cell|<frac|g|2>*<matrix|<tformat|<table|<row|<cell|0>>|<row|<cell|\h<rsup|2><rsub|i,r>+\h<rsup|2><rsub|i,f>-\h<rsup|2><rsub|i,r>-*\h<rsup|2><rsub|i,l>+\h<rsup|2><rsub|i,l>-\h<rsup|2><rsub|i,c>>>>>>=<frac|g|2>*<matrix|<tformat|<table|<row|<cell|0>>|<row|<cell|\h<rsup|2><rsub|i,f>*-\h<rsup|2><rsub|i,c>>>>>>>>>>
  </eqnarray*>

  It is thus tempting to rebalance the scheme using

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<cal-F\><rsub|l>>|<cell|=>|<cell|\<cal-F\><around|(|U<rsub|i+1/2->,U<rsub|i+1/2+>|)>+<frac|g|2>*<matrix|<tformat|<table|<row|<cell|0>>|<row|<cell|\h<rsup|2><rsub|i,r>-\h<rsup|2><rsub|i+1/2->+<around*|(|\h<rsub|i,r>+\h<rsub|i,f>|)>*<around*|(|z<rsub|i,r>-z<rsub|i>|)><with|color|dark
    green|-\h<rsup|2><rsub|i,f>>>>>>>,>>|<row|<cell|\<cal-F\><rsub|r>>|<cell|=>|<cell|\<cal-F\><around|(|U<rsub|i-1/2->,U<rsub|i-1/2+>|)>+<frac|g|2>*<matrix|<tformat|<table|<row|<cell|0>>|<row|<cell|\h<rsup|2><rsub|i,l>-\h<rsup|2><rsub|i-1/2+>+<around*|(|\h<rsub|i,l>+\h<rsub|i,c>|)>*<around*|(|z<rsub|i,l>-z<rsub|i>|)><with|color|dark
    green|-\h<rsup|2><rsub|i,c>>>>>>>,>>>>
  </eqnarray*>

  unfortunately this adds a term even in the case where <math|z=0>.
</body>

<\references>
  <\collection>
    <associate|audusse|<tuple|2|?>>
    <associate|auto-1|<tuple|1|?>>
    <associate|auto-2|<tuple|1.1|?>>
    <associate|auto-3|<tuple|1.2|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Balanced
      Saint-Venant scheme> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <with|par-left|<quote|1.5fn>|1.1<space|2spc>Uniform case
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>>

      <with|par-left|<quote|1.5fn>|1.2<space|2spc>Fine/coarse case
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3>>
    </associate>
  </collection>
</auxiliary>