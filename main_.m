%% plotting cheb
n=5;
tic
T=cheb(n);
toc
aantal_ev = 10000;
X = linspace(-1,1,aantal_ev);
figure()
hold on
for i=1:n+1
    plot(X,polyval(fliplr(T(i,:)),X))%Select the row in T, flip it to descending order and then evaluate in the given points X and plot.
end
ylim([-1 2]);
title('Chebychev veeltermen')
title(legend,'T_k')
legend(cellstr(num2str([0:n]', 'T_%-d')))

%% plotje van f_handle
f_handle = @(x) exp(x);
n = 6;
a = real(chebcoeff(f_handle,n));
a = a(1:round(length(a)/2));
T = cheb(n);
c = polyy(a,T);

aantal_ev = 100;
X = linspace(-1,1,aantal_ev);
y = zeros(1,aantal_ev);
for i=0:n
    y = y + c(i+1)*(X.^i);
end
    %ylabel('Amplitude');
    %xlabel('Frequentie (Hz)');


figure()
plot(X,y)
