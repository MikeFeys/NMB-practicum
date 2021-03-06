function [a] = chebcoeff(f_handle,N)
%evaluate op eenheidscirkel, punten van pi tot 2pi zijn identiek en dus
%nutteloos (na cos(theta)) daarom van 0 tot pi
evspace = cos(linspace(0,pi,N+1));

monsters = f_handle(evspace);

%uitbreiding (antiklokwijs, starten op x=1)     
monsters_eenheid_cirkel = [monsters,flip(monsters(2:end-1))];
        
%integraal met fft
fft_coeff = real(fft(monsters_eenheid_cirkel));
a = fft_coeff(1:N+1)/N; %tot N+1 maar anders hebben we er teveel

% laatste en eerste coeff delen door 2
a(1) = a(1)/2;
a(end) = a(end)/2;
end