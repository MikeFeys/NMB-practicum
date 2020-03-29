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
f_handle = @(x) x.^2;
n = 30;
a = chebcoeff(f_handle,n);
T = cheb(n);
c = poly(a,T);

aantal_ev = 100;
X = linspace(-1,1,aantal_ev);
y = zeros(1,aantal_ev);
for i=1:n+1
    y = y + c(i)*(X.^(i-1));
end
    %ylabel('Amplitude');
    %xlabel('Frequentie (Hz)')
figure()
plot(X,y)
if n==1
    title(['Chebychev veeltermbenadering van: ', strrep(char(f_handle),'@(x)','') ,' tot op 1ste orde.'])
else
    title(['Chebychev veeltermbenadering van: ', strrep(char(f_handle),'@(x)','') ,' tot op ',num2str(n),'de orde.'])
end
