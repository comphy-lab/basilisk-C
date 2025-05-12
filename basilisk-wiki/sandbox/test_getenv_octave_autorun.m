x=0:.1:10;
V = getenv('OCTAVE_AUTORUN')


if (V=='true')
plot(x,cos(x));
else
plot(x,sin(x));
end

print('-dsvg','-r80','test.svg');

%{
![this should be a cosine](test_getenv_octave_autorun/test.svg)
%}

