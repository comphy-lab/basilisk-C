<TeXmacs|1.0.7.15>

<style|generic>

<\body>
  We suppose that the number of grid points necessary to describe a given
  phenomenon scales like

  <\equation*>
    \<Delta\><rsup|-d>
  </equation*>

  with <math|\<Delta\>> the grid size and <math|d> the fractal dimension of
  the process. By construction

  <\equation*>
    \<Delta\><rsub|l>\<sim\>2<rsup|-l>
  </equation*>

  with <math|l> the level of the quadtree cell. The number <math|n<rsub|l>>
  of grid points at level <math|l> thus verifies

  <\equation*>
    2<rsup|d*l>=2<rsup|d<around*|(|l-1|)>>-<frac|n<rsub|l>|2<rsup|D>>+n<rsub|l>
  </equation*>

  with <math|D> the dimension of the containing space. This gives

  <\equation*>
    n<rsub|l>=<frac|2<rsup|d*l>*<around*|(|1-2<rsup|-d>|)>|1-<frac|1|2<rsup|D>>>=<frac|2<rsup|D>|2<rsup|D>-1>*2<rsup|d*l>*<around*|(|1-2<rsup|-d>|)>
  </equation*>

  and the fractal dimension verifies

  <\equation*>
    d=log<rsub|2><around*|(|<frac|n<rsub|l+1>|n<rsub|l>>|)>
  </equation*>

  If a constant timestep <math|\<Delta\>t\<sim\>\<Delta\><rsub|L>>, with
  <math|L> the maximum level, is chosen to integrate the entire system, the
  cost will scale like <math|2<rsup|L<around*|(|d+1|)>>>. If a multirate time
  integration is used, the integration cost per level scales like
  <math|2<rsup|l>*n<rsub|l>=2<rsup|l<around*|(|d+1|)>>*<around*|(|1-2<rsup|-d>|)>>
  and the total integration cost scales like

  <\equation*>
    <big|sum><rsub|l>2<rsup|l>*n<rsub|l>=<frac|2<rsup|D>|2<rsup|D>-1>*<around*|(|1-2<rsup|-d>|)>*<big|sum><rsub|l><around*|(|2<rsup|d+1>|)><rsup|l>==<frac|2<rsup|D>|2<rsup|D>-1>*<frac|1-2<rsup|-d>|2<rsup|d+1>-1>*<around*|(|2<rsup|<around*|(|L+1|)>*<around*|(|d+1|)>>-1|)>
  </equation*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|<big|sum><rsub|l>2<rsup|l>*n<rsub|l>>|<cell|=>|<cell|<frac|2<rsup|D>|2<rsup|D>-1>*<around*|(|1-2<rsup|-d>|)>*<big|sum><rsub|l><around*|(|2<rsup|d+1>|)><rsup|l>>>|<row|<cell|>|<cell|=>|<cell|<frac|2<rsup|D>|2<rsup|D>-1>*<around*|(|1-2<rsup|-d>|)>*<frac|1-2<rsup|<around*|(|L+1|)>*<around*|(|d+1|)>>|1-2<rsup|d+1>>>>|<row|<cell|>|<cell|=>|<cell|<frac|2<rsup|D>|2<rsup|D>-1>*<frac|1-2<rsup|-d>|2<rsup|d+1>-1>*<around*|(|2<rsup|<around*|(|L+1|)>*<around*|(|d+1|)>>-1|)>>>>>
  </eqnarray*>

  The efficiency gain, relative to using a constant timestep is thus

  <\equation*>
    <frac|2<rsup|L<around*|(|d+1|)>>|<big|sum><rsub|l>2<rsup|l>*n<rsub|l>>=<frac|2<rsup|D>-1|2<rsup|D>>*<frac|2<rsup|L<around*|(|d+1|)>>*<around*|(|2<rsup|d+1>-1|)>|<around*|(|1-2<rsup|-d>|)>*<around*|(|2<rsup|<around*|(|L+1|)>*<around*|(|d+1|)>>-1|)>>
  </equation*>

  If <math|2<rsup|<around*|(|L+1|)>*<around*|(|d+1|)>>\<gg\>1>, this gives

  <\equation*>
    <frac|2<rsup|L<around*|(|d+1|)>>|<big|sum><rsub|l>2<rsup|l>*n<rsub|l>>\<simeq\><frac|2<rsup|D>-1|2<rsup|D>>*<frac|2<rsup|d+1>-1|2<rsup|d+1>-2<rsup|*>>
  </equation*>

  Note that this gain depends only on the fractal dimension, not on the
  maximum level <math|L>.

  <big-figure|<image|multirate.eps|0.7par|||>|Multirate efficiency gain as a
  function of fractal dimension for <math|D=2>.>
</body>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|figure>
      <tuple|normal|Multirate efficiency gain as a function of fractal
      dimension.|<pageref|auto-1>>
    </associate>
  </collection>
</auxiliary>