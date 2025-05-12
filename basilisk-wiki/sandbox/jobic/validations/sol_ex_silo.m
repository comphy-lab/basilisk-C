function [x] = sol_ex_silo(pval,posx,pmilieu)
  dg=0.04;
  mus=0.383;
  dmu=0.26;
  I0=0.279;
  x0=mus*pval;
  dmuP=dmu*pval;
  x=zeros(1,length(posx));

  for i=1:length(posx)
    xx=abs(posx(i));
    if (xx <  x0) 
      x(i)=-I0*sqrt(pval)*((x0-pmilieu)+dmuP*log(dmuP/(dmuP-(pmilieu-x0))));
      ((x0-pmilieu)+dmuP*log(dmuP/(dmuP-(pmilieu-x0))))
      I0*sqrt(pval)
    else 
      x(i)=-I0*sqrt(pval)*((xx-pmilieu)+dmuP*log((dmuP-(xx-x0))/(dmuP-(pmilieu-x0))));
    end
  end


end
