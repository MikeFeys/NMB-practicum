%% 1.2.3 f(x) benaderen + convergentiesnelheid enzo
clear
tic
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
TijdTerInfoInSeconden=toc
%% plotting cheb
clear

tic
n=4;
aantal_ev = 100;

T=cheb(n)
X = linspace(-1,1,aantal_ev);
figure()
hold on
set(gca,'colororder',[0 0 1; 0 1 0],'linestyleorder',{'-','-.',':','--','-*',':s','--^'},'nextplot','add')
for i=1:n+1
    plot(X,polyval(fliplr(T(i,:)),X))%Select the row in T, flip it to descending order and then evaluate in the given points X and plot.
end
ylim([-1 2]);
title(['Chebyshev veeltermen tot ',num2str(n) , ' orde'])
title(legend,'T_k')
legend(cellstr(num2str([0:n]', 'T_%-d')))
TijdTerInfoInSeconden=toc
%% plotje van f_handle
clear
tic
%f_handle = @(x) (x-1)./(1+6*x.^2);
%f_handle = @(x) log(x+2).*sin(10*x);
f_handle = @(x) x.^5-x.^4+x.^3-x.^2+x-1;

aantal_ev = 100;
X = linspace(-1,1,aantal_ev);
figure()
hold on
fplot(f_handle,[-1 1])
set(gca,'colororder',[0 0 1; 0 1 0],'linestyleorder',{'-','-.',':','--','-*',':s','--^'},'nextplot','add')
for n =[2 4 9 11 30] %Aantal verschillende benaderingen plotten
a = chebcoeff(f_handle,n);
T = cheb(n);
c = poly(a,T);
y = zeros(1,aantal_ev);
for i=1:n+1
    y = y + c(i)*(X.^(i-1));
end
plot(X,y,'DisplayName',['Benadering van ',num2str(n),' orde'])
hold on
legend()
title(['Chebychev veeltermbenadering van: ', strrep(char(f_handle),'@(x)','') ,' tot op ',num2str(n),' orde.'])
end
TijdTerInfoInSeconden=toc
%% Convergentiesnelheid LUKAS
clear
tic
f_handle = @(x) (x-1)./(1+6*x.^2);
%f_handle = @(x) log(x+2).*sin(10*x);
%f_handle = @(x) x.^5-x.^4+x.^3-x.^2+x-1;
max_ord = 60;

X = linspace(-1,1,100);
max_fout = zeros(1,max_ord);
T = cheb(max_ord);
figure()
for n = 1:max_ord
    
a = chebcoeff(f_handle,n);
c = poly(a,T(1:n+1,1:n+1));
syms x;
y(x) = 0.*x;

for i=1:n+1
    y = y + c(i)*(x.^(i-1));
end

[x_min,fval] = fminsearch(-abs(f_handle-y),0);
    if x_min>1 || x_min<-1
        disp("ERROR: Minima is outside boundaries of [-1,1]");
    end
max_fout(n) = -fval;
%hold on
%fplot(abs(f_handle-y),[-1 1],'DisplayName',['Fout van orde ',num2str(n)])
end
%fplot(log(abs(f_handle-y)),[-1 1])
%legend()
hold off
%figure()
plot(1:n,max_fout)
TijdTerInfoInSeconden=toc
%% BATS
clear
tic
f_handle = @(x) (x-1)./(1+6*x.^2);
%f_handle = @(x) log(x+2).*sin(10*x);
%f_handle = @(x) x.^5-x.^4+x.^3-x.^2+x-1;
max_ord = 60;

max_fout = zeros(1,max_ord);
T = cheb(max_ord);
for n = 1:max_ord
    
a = chebcoeff(f_handle,n);
c = poly(a,T(1:n+1,1:n+1));
y = poly2sym(fliplr(c));

%limatation of fminbnd is that it may be stuck in local minima so evaluate
%max error per interval and seek maxiumum error of all the intervals
evspace = fliplr(cos(linspace(0,pi,n+1)) );
max_intval = zeros(0,length(evspace)-1); 
for i=2:length(evspace)
    [~,fval,~,~] = fminbnd( matlabFunction(-abs(f_handle-y)),evspace(i-1),evspace(i));
    max_intval(i) = -fval;
end
max_fout(n) = max(max_intval);
end
figure()
plot(1:n,max_fout)

TijdTerInfoInSeconden=toc
%%
%Gemiddelde fout convergentie (mike)
clear
tic
%f_handle = @(x) (x-1)./(1+6*x.^2);
%f_handle = @(x) log(x+2).*sin(10*x);
f_handle = @(x) x.^5-x.^4+x.^3-x.^2+x-1;

max_orde=65;%Maximum chebyshev orde benadering
aantal_ev = 100;%Aantal interpolatiepunten

X = linspace(-1,1,aantal_ev);%Gelijk verdeelde interpolatiepunten maken
Gemiddeldefoutlog=zeros(1,max_orde);%Preallocate voor snelheid
Functiewaarden=arrayfun(f_handle,X);%Alle functiewaarden van de effectieve functie over de interpolatiepunten
T = cheb(max_orde);%Geeft alle veeltermen in matrix tot max_orde dus kan op voorhand al berekend worden en dan slechts opvragen

for n =1:max_orde 
a = chebcoeff(f_handle,n);
c = poly(a,T(1:n+1,1:n+1));
y = zeros(1,aantal_ev);

for i=1:n+1
    y = y + c(i)*(X.^(i-1));
end
%Ik neem aan dat sum al alles ordend en klein naar groot de som neemt om de
%kleineste fout te krijgen
Gemiddeldefoutlog(n)=log10(sum(abs(Functiewaarden-y))/aantal_ev);%Log van de som van de absolute waarden van het verschil van effectieve waarde en benaderingswaarde gedeeld door het aantal interpolatiepunten 
end
plot(1:max_orde,Gemiddeldefoutlog,'HandleVisibility','off')
hold on
plot(1:max_orde,Gemiddeldefoutlog,'b*')
xlabel('Orde')
ylabel('log(gemiddelde fout)')
legend('gemiddelde fout per orde')
title(['Convergentie van ', strrep(char(f_handle),'@(x)','')])
TijdTerInfoInSeconden=toc