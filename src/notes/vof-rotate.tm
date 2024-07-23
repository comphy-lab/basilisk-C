<TeXmacs|1.99.14>

<style|generic>

<\body>
  <\big-figure|<with|gr-mode|<tuple|group-edit|zoom>|gr-frame|<tuple|scale|0.840896cm|<tuple|0.5gw|0.5gh>>|gr-geometry|<tuple|geometry|0.5par|0.6par|center>|gr-grid|<tuple|empty>|gr-grid-old|<tuple|cartesian|<point|0|0>|1>|gr-edit-grid-aspect|<tuple|<tuple|axes|none>|<tuple|1|none>|<tuple|10|none>>|gr-edit-grid|<tuple|empty>|gr-edit-grid-old|<tuple|cartesian|<point|0|0>|1>|magnify|0.840896415|gr-arrow-end|\<gtr\>|gr-arrow-begin|\<less\>|gr-auto-crop|true|<graphics||<line|<point|-2|2>|<point|2.0|2.0>|<point|2.0|-2.0>|<point|-2.0|-2.0>|<point|-2.0|2.0>>|<with|dash-style|11100|arrow-end|\<gtr\>|arrow-begin|\<less\>|<line|<point|0|3>|<point|0.0|0.0>|<point|3.0|0.0>>>|<text-at|<math|y>|<point|0.2|3>>|<text-at|<math|x>|<point|3.0|-0.4>>|<with|line-width|2ln|<line|<point|-2.4|3.0>|<point|1.8|-3.0>>>|<with|dash-style|10|line-width|2ln|<line|<point|-1.0|3.0>|<point|0.4|-3.0>>>|<text-at|<math|<frac|1|2>>|<point|2.2|-0.6>>|<text-at|<math|-<frac|1|2>>|<point|-2.7|-0.6>>|<text-at|<math|<frac|1|2>>|<point|0.2|2.3>>|<text-at|<math|-<frac|1|2>>|<point|-0.7|-2.6>>|<text-at|<math|O>|<point|-0.7|-0.3>>|<text-at|<math|x<rsub|1>>|<point|-1.9|1.6>>|<text-at|<math|x<rsub|2>>|<point|1.1|-1.9>>|<with|arrow-end|\<gtr\>|arrow-begin|\<less\>|<line|<point|-1.7|2.1>|<point|-0.900000000000001|2.11366612997792>>>|<with|arrow-end|\<gtr\>|arrow-begin|\<less\>|<line|<point|0.3|-2.1>|<point|1.1|-2.1>>>|<text-at|<math|r>|<point|-1.4|2.2>>|<text-at|<math|r>|<point|0.6|-2.4>>|<text-at||<point|2.7|-1.3>>>>>
    \;
  </big-figure>

  The equation for the interface (black line) is

  <\equation*>
    n<rsub|x>*x+n<rsub|y>*y=\<alpha\>
  </equation*>

  Using the fact that both interfaces intersect in <math|O> on the
  <math|x>-axis, the equation for the rotated interface (dashed line) can be
  written

  <\equation*>
    n<rsub|x>*x+n<rsup|1><rsub|y>*y=\<alpha\>
  </equation*>

  The interfaces verify

  <\eqnarray*>
    <tformat|<table|<row|<cell|n<rsub|x>*x<rsub|1>+<frac|n<rsub|y>|2>>|<cell|=>|<cell|\<alpha\>>>|<row|<cell|n<rsub|x>*x<rsub|2>-<frac|n<rsub|y>|2>>|<cell|=>|<cell|\<alpha\>>>|<row|<cell|n<rsub|x>*<around*|(|x<rsub|1>+r|)>+<frac|n<rsup|1><rsub|y>|2>>|<cell|=>|<cell|\<alpha\>>>|<row|<cell|n<rsub|x>*<around*|(|x<rsub|2>-r|)>-<frac|n<rsup|1><rsub|y>|2>>|<cell|=>|<cell|\<alpha\>>>>>
  </eqnarray*>

  which gives

  <\eqnarray*>
    <tformat|<table|<row|<cell|n<rsup|1><rsub|y>>|<cell|=>|<cell|n<rsub|y>-2*n<rsub|x>*r>>>>
  </eqnarray*>

  The displacement due to rotation <math|r> is given by the horizontal shear
  and can be written

  <\equation*>
    r=<frac|\<Delta\>t|2>*<frac|\<partial\>u|\<partial\>y>
  </equation*>

  so that the equation of the rotated interface is simply

  <\equation*>
    n<rsub|x>*x+<around*|(|n<rsub|y>-n<rsub|x>*\<Delta\>t*<frac|\<partial\>u|\<partial\>y>|)>*y=\<alpha\>
  </equation*>
</body>

<\initial>
  <\collection>
    <associate|page-medium|paper>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|1|../../../.TeXmacs/texts/scratch/no_name_1.tm>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|figure>
      <tuple|normal|<\surround|<hidden-binding|<tuple>|1>|>
        \;
      </surround>|<pageref|auto-1>>
    </associate>
  </collection>
</auxiliary>