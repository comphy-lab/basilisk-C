fs = fopen ("fft_static.dat", "w");

for Re = 250:500:750

    # Load file	
    fp = fopen (strcat("../../strouhal_rectangle/force-Re",num2str(Re)), "r");
    data = dlmread (fp,' ', 0, 1);

    # Extract data and remove transient (before t=80)
    tt = data (:,1);

    it = 1;
    while (tt(it) < 80)
      it++;
   endwhile

   t = data (it:end,1);
   dt = data (it:end,2);
   Cd = data (it:end,3) + data (it:end,5);
   Cl = data (it:end,4) + data (it:end,6);

   # Resample the data with uniform samplying frequency
   N = length(t); 	    	 	   	   % Signal Length
   tr = linspace(min(t), max(t), N);  	   	   % Uniformly-Sampled Time Vector
   Clr = interp1(t,Cl,tr,'spline');       	   % Resampled Signal Vector

   # Apply fft to Cl
   FT = fft(Clr,N);				   % Fourier Transform
   FTS = abs(fftshift(FT));			   % Shifted FFT

   # Define the frequency domain
   Ts = mean(diff(tr));				   % Sampling Interval
   Fs = 1/Ts;					   % Sampling Frequency
   Fn = Fs/2;					   % Nyquist Frequency
   dF = Fs/N;
   F = [-Fn:dF:Fn - dF];

   # Sort the signal
   [M, I] = sort (FTS,"descend");

   # Compute the weighted average of the frequency
   tresh = 55.
   Freq = 0.;
   sum_M = 0.;
   for i=1:2
       if (M(i) > tresh)
       	  F(I(i)), M(i)	
       	  Freq += abs(F(I(i)))*M(i);
       	  sum_M += M(i);
       endif
   endfor
   Freq /= sum_M;

   # Print the data into a file
   fprintf(fs, "%g %g\n", Re, Freq);
   fclose (fp);
   
endfor
fclose (fs);