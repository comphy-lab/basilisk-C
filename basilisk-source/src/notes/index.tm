<TeXmacs|2.1.1>

<style|<tuple|generic|maxima>>

<\body>
  <section*|Indexing for <verbatim|grid/multigrid.h>>

  <subsection|1D>

  Given <math|l> the grid level, the number of grid points along one
  dimension is

  <\equation*>
    s<around*|(|l|)>=2<rsup|l>+2*g
  </equation*>

  with <math|g> the number of ghost layers (usually 2).

  The index of field <math|a> at coordinates <math|<around*|(|i|)>> is given
  by

  <\equation*>
    i+a*s<around*|(|l|)>
  </equation*>

  In a similar fashion, the index of field <math|a> at coordinates
  <math|<around*|(|i|)>> and level <math|l> is given by

  <\equation*>
    i+l*<big|sum><rsup|k\<less\>l><rsub|k=0>s<around*|(|k|)>+a*<big|sum><rsup|k\<less\>depth+1><rsub|k=0>s<around*|(|k|)>
  </equation*>

  We can further develop

  <\equation*>
    <big|sum><rsup|k\<less\>l><rsub|k=0>s<around*|(|k|)>=<big|sum><rsup|k\<less\>l><rsub|k=0><around*|(|2<rsup|k>+2*g|)>
  </equation*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|<big|sum><rsup|k\<less\>l><rsub|k=0>s<around*|(|k|)>>|<cell|=>|<cell|<big|sum><rsup|k\<less\>l><rsub|k=0><around*|(|2<rsup|k>+2*g|)>>>|<row|<cell|>|<cell|=>|<cell|2<rsup|l>-1+2*g*l>>>>
  </eqnarray*>

  <subsection|2D>

  Given <math|l> the grid level, the number of grid points along one
  dimension is

  <\equation*>
    s<around*|(|l|)>=2<rsup|l>+2*g
  </equation*>

  with <math|g> the number of ghost layers (usually 2).

  If we choose to pack elements row-by-row (i.e. elements are closest along
  rows), the index of field <math|a> at coordinates <math|<around*|(|i,j|)>>
  is given by

  <\equation*>
    j+i*s<around*|(|l|)>+a*s<around*|(|l|)><rsup|2>
  </equation*>

  In a similar fashion, if we pack row-by-row and level-by-level, the index
  of field <math|a> at coordinates <math|<around*|(|i,j|)>> and level
  <math|l> is given by

  <\equation*>
    j+i*s<around*|(|l|)>+l*<big|sum><rsup|k\<less\>l><rsub|k=0>s<around*|(|k|)><rsup|2>+a*<big|sum><rsup|k\<less\>depth+1><rsub|k=0>s<around*|(|k|)><rsup|2>
  </equation*>

  We can further develop

  <\equation*>
    <big|sum><rsup|k\<less\>l><rsub|k=0>s<around*|(|k|)><rsup|2>=<big|sum><rsup|k\<less\>l><rsub|k=0><around*|(|2<rsup|k>+2*g|)><rsup|2>
  </equation*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|<big|sum><rsup|k\<less\>l><rsub|k=0>s<around*|(|k|)><rsup|2>>|<cell|=>|<cell|<big|sum><rsup|k\<less\>l><rsub|k=0><around*|(|2<rsup|k>+2*g|)><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsup|k\<less\>l><rsub|k=0>4<rsup|*k>+4*g*2<rsup|k>+4*g<rsup|2>>>|<row|<cell|>|<cell|=>|<cell|<frac|4<rsup|l>-1|3>+4*g*<around*|(|2<rsup|l>-1+g*l|)>>>>>
  </eqnarray*>

  <subsection|3D>

  <\equation*>
    k+s<around*|(|l|)>*<around*|(|j*+i*s<around*|(|l|)>|)>+l*<big|sum><rsup|k\<less\>l><rsub|k=0>s<around*|(|k|)><rsup|3>+a*<big|sum><rsup|k\<less\>depth+1><rsub|k=0>s<around*|(|k|)><rsup|3>
  </equation*>

  We can further develop

  <\eqnarray*>
    <tformat|<table|<row|<cell|<big|sum><rsup|k\<less\>l><rsub|k=0>s<around*|(|k|)><rsup|3>>|<cell|=>|<cell|<big|sum><rsup|k\<less\>l><rsub|k=0><around*|(|2<rsup|k>+2*g|)><rsup|3>>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsup|k\<less\>l><rsub|k=0>3*g*2<rsup|2*k+1>+3*g<rsup|2>*2<rsup|k+2>+2<rsup|3*k>+8*g<rsup|3>>>|<row|<cell|>|<cell|=>|<cell|<frac|8<rsup|l>-1|7>+2*g*<space|0.17em><around*|(|4<rsup|l>-1|)>+12*g<rsup|2>*<around*|(|2<rsup|l>-1|)>+8*g<rsup|3>*l>>>>
  </eqnarray*>
</body>

<\initial>
  <\collection>
    <associate|page-height|auto>
    <associate|page-type|letter>
    <associate|page-width|auto>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|?|?>>
    <associate|auto-2|<tuple|1|?>>
    <associate|auto-3|<tuple|2|?>>
    <associate|auto-4|<tuple|3|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Indexing
      for <with|font-family|<quote|tt>|language|<quote|verbatim>|grid/multigrid.h>>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <with|par-left|<quote|1tab>|1<space|2spc>1D
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>>

      <with|par-left|<quote|1tab>|2<space|2spc>2D
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3>>

      <with|par-left|<quote|1tab>|3<space|2spc>3D
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4>>
    </associate>
  </collection>
</auxiliary>