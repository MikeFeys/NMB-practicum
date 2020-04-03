%% 1.2.3 f(x) benaderen + convergentiesnelheid enzo
syms x;
f1(x) = (x-1)./(1+6.*x.^2);

f_handle = @(x) (x-1)./(1+6.*x.^2);
n_max =100;
T = cheb(n_max);

aantal_ev = 100;
X = linspace(-1,1,aantal_ev);
for n=2:n_max
    a = chebcoeff(f_handle,n);
    c = poly(a,T);
    %polyval(fliplr(c),X)
    
    %plot(X, abs((X-1)./(1+6.*X.^2)-polyval(fliplr(c),X)));
end
    plot(X, abs((X-1)./(1+6.*X.^2)-polyval(fliplr(c),X)));



%%
% plotten
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
f_handle = @(x) (x-1)/(1+6*x.^2);

aantal_ev = 100;
X = linspace(-1,1,aantal_ev);
figure()
hold on
fplot(f_handle,[-1 1])
set(gca,'colororder',[0 0 1; 0 1 0],'linestyleorder',{'-','-.',':','--','-*',':s','--^'},'nextplot','add')
for n =[2 4 6 8 12 16 20] %Aantal verschillende benaderingen plotten
a = chebcoeff(f_handle,n);
T = cheb(n);
c = poly(a,T);
y = zeros(1,aantal_ev);
for i=1:n+1
    y = y + c(i)*(X.^(i-1));
end
 if n==1
    ordenaam=' ste orde' ;
else 
    ordenaam=' de orde' ;
 end
plot(X,y,'DisplayName',['Benadering van ',num2str(n),ordenaam])
hold on
legend()
if n==1
    title(['Chebychev veeltermbenadering van: ', strrep(char(f_handle),'@(x)','') ,' tot op 1ste orde.'])
else
    title(['Chebychev veeltermbenadering van: ', strrep(char(f_handle),'@(x)','') ,' tot op ',num2str(n),'de orde.'])
end


end

%% Convergentiesnelheid
max_ord = 100;
max_fout = zeros(1,max_ord);

f_handle = @(x) (x-1)/(1+6*x.^2);

 

for n = 1:max_ord
    
a = chebcoeff(f_handle,n);
T = cheb(n);
c = poly(a,T);
syms x;
y(x) = 0.*x;

for i=1:n+1
    y = y + c(i)*(x.^(i-1));
end

max_fout(1,n) = max();
end
figure()
plot(1:n,max_fout)









