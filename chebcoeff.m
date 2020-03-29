function [a] = chebcoeff(f_handle,N)
%evaluate
evspace = linspace(0,1,N+1);
monsters = f_handle(cos(pi.*evspace./(N+1)));

%uitbreiding
flipmonsters = flip(monsters(2:end-1));
monsters = [monsters flipmonsters];

%integraal
a = abs(fft(monsters)./(N+1));
end