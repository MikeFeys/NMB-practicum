%% Plotten van enkele Chebyshev-veeltermen om te testen
clear
tic

%Invoer parameters.
n=4;
aantal_ev = 100;

%Bepalen interpolatiepunten en Chebyshev-veeltermen opvragen.
T=cheb(n);
X = linspace(-1,1,aantal_ev);

%Plot aanpassen.
figure()
hold on
set(gca,'colororder',[0 0 1; 0 1 0],'linestyleorder',{'-','-.',':','--','-*',':s','--^'},'nextplot','add')

%De rij in T in dalende orde zetten en evalueren in de interpolatiepunten en plotten.
for i=1:n+1
    plot(X,polyval(fliplr(T(i,:)),X))
end

%Plot aanpassen.
ylim([-1 2]);
title(['Chebyshev veeltermen tot ',num2str(n) , ' orde'])
title(legend,'T_k')
legend(cellstr(num2str([0:n]', 'T_%-d')))

%Ter info de tijd die de berekening nam.  
TijdTerInfoInSeconden=toc;
disp(['Het uitvoeren van alle code duurde: ', num2str(TijdTerInfoInSeconden),' seconden.'])
%% Plotten van functies met verschillende benaderingen
clear
tic
%Invoeren parameters.
aantal_ev = 100;
LijstBenaderingsOrdes=[2 4 9 11 30];

%Afhankelijk van welke functie de rest in commentaar zetten.
%f_handle = @(x) (x-1)./(1+6*x.^2);
%f_handle = @(x) log(x+2).*sin(10*x);
f_handle = @(x) x.^5-x.^4+x.^3-x.^2+x-1;

%Interpolatiepunten bepalen.
X = linspace(-1,1,aantal_ev);

%Plot aanpassen.
figure()
hold on
fplot(f_handle,[-1 1])
set(gca,'colororder',[0 0 1; 0 1 0],'linestyleorder',{'-','-.',':','--','-*',':s','--^'},'nextplot','add')

%Matrix chebychev veeltermen van maximale orde op voorhand berekenen en daaruit voor kleinere orde selecteren.
T = cheb(max(LijstBenaderingsOrdes));

%Voor de ordes in LijstBenaderingsOrdes de orde plotten.
for n =LijstBenaderingsOrdes 
    
%Opvragen waardes.
Coefficientenak = chebcoeff(f_handle,n);
Coefficientenyn = poly(Coefficientenak,T(1:n+1,1:n+1));

%Preallocate voor snelheid.
y = zeros(1,aantal_ev);

%y-waarden berekenen.
for i=1:n+1
    y = y + Coefficientenyn(i)*(X.^(i-1));
end

%Plot.
plot(X,y,'DisplayName',['Benadering van ',num2str(n),' orde'])
hold on
legend()
title(['Chebychev veeltermbenadering van: ', strrep(char(f_handle),'@(x)','') ,' tot op ',num2str(n),' orde.'])
end

%Ter info de tijd die de berekening nam.  
TijdTerInfoInSeconden=toc;
disp(['Het uitvoeren van alle code duurde: ', num2str(TijdTerInfoInSeconden),' seconden.'])
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
%Ter info de tijd die de berekening nam.  
TijdTerInfoInSeconden=toc;
disp(['Het uitvoeren van alle code duurde: ', num2str(TijdTerInfoInSeconden),' seconden.'])
%% BATS
clear
tic
%f_handle = @(x) (x-1)./(1+6*x.^2);
f_handle = @(x) log(x+2).*sin(10*x);
%f_handle = @(x) x.^5-x.^4+x.^3-x.^2+x-1;
max_ord = 100;

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
plot(1:n,log10(max_fout),'HandleVisibility','off')
hold on
plot(1:max_ord,log10(max_fout),'b*')
xlabel('Orde')
ylabel('log(maximale fout)')
legend('maximale fout per orde')
title(['Convergentie van ', strrep(char(f_handle),'@(x)','')])

%Ter info de tijd die de berekening nam.  
TijdTerInfoInSeconden=toc;
disp(['Het uitvoeren van alle code duurde: ', num2str(TijdTerInfoInSeconden),' seconden.'])
%% Convergentie via gemiddelde fout (Mike)
clear
tic

%Afhankelijk van welke functie de andere 2 in commentaar zetten.
f_handle = @(x) (x-1)./(1+6*x.^2);
%f_handle = @(x) log(x+2).*sin(10*x);
%f_handle = @(x) x.^5-x.^4+x.^3-x.^2+x-1;

%Parameters ingeven.
max_orde=65;%Maximum chebyshev orde benadering
aantal_ev = 100;

%Interpolatiepunten bepalen.
X = linspace(-1,1,aantal_ev);
%Preallocate voor snelheid
Gemiddeldefoutlog=zeros(1,max_orde);
%Alle functiewaarden van de effectieve functie over de interpolatiepunten
Functiewaarden=arrayfun(f_handle,X);
%Matrix chebychev veeltermen van maximale orde op voorhand berekenen en daaruit voor kleinere orde selecteren.
T = cheb(max_orde);

for n =1:max_orde 

%Opvragen waardes.
Coefficientenak = chebcoeff(f_handle,n);
Coefficientenyn = poly(Coefficientenak,T(1:n+1,1:n+1));

%Preallocate voor snelheid.
y = zeros(1,aantal_ev);

%y-waarden berekenen.
for i=1:n+1
    y = y + Coefficientenyn(i)*(X.^(i-1));
end
%Ik neem aan dat sum de som van klein naar groot neemt om de kleinste fout
%te krijgen.
%Log van de som van de absolute waarden van het verschil van effectieve waarde
%en benaderingswaarde gedeeld door het aantal interpolatiepunten om het
%log(gemiddeld) te krijgen
Gemiddeldefoutlog(n)=log10(sum(abs(Functiewaarden-y))/aantal_ev); 
end

%Plot.
plot(1:max_orde,Gemiddeldefoutlog,'HandleVisibility','off')
hold on
plot(1:max_orde,Gemiddeldefoutlog,'b*')
xlabel('Orde')
ylabel('log10(gemiddelde fout)')
legend('gemiddelde fout per orde')
title(['Convergentie van ', strrep(char(f_handle),'@(x)','')])

%Ter info de tijd die de berekening nam.  
TijdTerInfoInSeconden=toc;
disp(['Het uitvoeren van alle code duurde: ', num2str(TijdTerInfoInSeconden),' seconden.'])
%% Deel 2: Bivariate kleinste-kwadratenveeltermbenadering
%nuttige Matlab-functies zijn meshgrid; polyval2; scatter3 en surf
clear
tic
%Parameters ingeven.
m=7;
n=7;
M=31;
N=31;
f_handle= @(x1,y1) sin((2*x1-1).^2+2*y1);

%Bepalen vectoren waarvan het cartesisch product puntenrooster vormt.
if M==N
    x=linspace(-1,1,M);
    y=x;
else
    x=linspace(-1,1,M);
    y=linspace(-1,1,N);
end
%Puntenrooster maken
[X,Y]=meshgrid(x,y);
%Functiewaarden van gegeven functie op puntenrooster.
F=f_handle(X,Y);

%Plot.
figure()
surf(X,Y,F)
xlabel('X')
ylabel('Y')
title('Werkelijke 3D plot van sin((2x-1)^2+2y)')

%Voor f(x,y)=sin((2*x-1)^2+2y)
    
    %Waarden opvragen.
    CoefMatrix=kkb2d(x, y, F, m, n);
    z=polyval2(CoefMatrix,x,y);%Approximatie.
    
    %Plot.
    figure()
    surf(X,Y,z,'HandleVisibility','off')
    hold on
    scatter3(reshape(X,[1,M*N]),reshape(Y,[1,M*N]),reshape(z,[1,M*N]),'filled','black')
    xlabel('X')
    ylabel('Y')
    legend('Datapunten')
    title(['Benaderende 3D plot van sin((2x-1)^2+2y) tot op graad m=',num2str(m),' en n=',num2str(n)])
    
%Voor F=membrane(1)
    %Plot.
    figure()
    surf(X,Y,membrane(1))
    xlabel('X')
    ylabel('Y')
    title('Werkelijke 3D plot van membrane(1)')
    
    %Waarden opvragen.
    CoefMatrix=kkb2d(x, y, membrane(1), m, n);
    z=polyval2(CoefMatrix,x,y);%Approximatie.
    
    %Plot.
    figure()
    surf(X,Y,z,'HandleVisibility','off')
    hold on
    scatter3(reshape(X,[1,M*N]),reshape(Y,[1,M*N]),reshape(z,[1,M*N]),'filled','black')
    xlabel('X')
    ylabel('Y')
    legend('Datapunten')
    title(['Benaderende 3D plot van functie membrane(1) tot op graad m=',num2str(m),' en n=',num2str(n)])
    
%Ter info de tijd die de berekening nam.  
TijdTerInfoInSeconden=toc;
disp(['Het uitvoeren van alle code duurde: ', num2str(TijdTerInfoInSeconden),' seconden.'])
%% Benaderingsfout plotten
clear
tic
%Parameters ingeven.
M=31;
N=31;
f_handle= @(x1,y1) sin((2*x1-1).^2+2*y1);

%Bepalen vectoren waarvan het cartesisch product puntenrooster vormt.
x=linspace(-1,1,M);
y=linspace(-1,1,N);

%Puntenrooster maken.
[X,Y]=meshgrid(x,y);
%Functiewaarden van gegeven functie op puntenrooster.
F=f_handle(X,Y);
%Preallocate voor snelheid.
Foutenlijst=zeros(2,20);

%Voor beide functies meteen de fouten berekenen.
for m=1:20
n=m;
CoefMatrix1=kkb2d(x, y, F, m, n);
CoefMatrix2=kkb2d(x, y, membrane(1), m, n);

z1=polyval2(CoefMatrix1,x,y);
z2=polyval2(CoefMatrix2,x,y);

Foutenlijst(1,m)=norm(F-z1);
Foutenlijst(2,m)=norm(F-z2);
end

%Plot.
plot(1:20,Foutenlijst(1,:))
hold on
plot(1:20,Foutenlijst(2,:))
xlabel('Graad van benadering')
ylabel('Norm van de benaderingsfout')
legend('Benaderingsfout van sin((2x-1)^2+2y)','Benaderingsfout van de functie membrane(1)')
title('Plot van de benaderingsfouten tot op orde m=20 en n=20')

%Ter info de tijd die de berekening nam.  
TijdTerInfoInSeconden=toc;
disp(['Het uitvoeren van alle code duurde: ', num2str(TijdTerInfoInSeconden),' seconden.'])
%% Benadering functie van Runge
clear
tic

%Parameters ingeven.
M=30;
N=M;
m=7;% This needs to be 7??
n=7;
f_handle= @(x1,y1) 1./(1+25*x1.^2+25*y1.^2);

%Lijst voor maximale afwijking
MaxFoutBenadering=zeros(2,1);
%Bepalen vectoren waarvan het cartesisch product puntenrooster vormt.
x=linspace(-1,1,M);
y=x;
%Puntenrooster maken.
[X,Y]=meshgrid(x,y);
%Functiewaarden van gegeven functie op puntenrooster.
F=f_handle(X,Y);

%Plot.
figure()
surf(X,Y,F)
xlabel('X')
ylabel('Y')
title('Werkelijke 3D plot van functie van Runge met equidistante verdeling')

%Equidistante verdeling

    %Waarden opvragen.
    CoefMatrix=kkb2d(x, y, F, m, n);
    z=polyval2(CoefMatrix,x,y);%Approximatie
    %Maximale fout bepalen.
    MaxFoutBenadering(1,1)=max(abs(F-z),[], 'all');
    
    %Plot.
    figure()
    surf(X,Y,z,'HandleVisibility','off')
    hold on
    scatter3(reshape(X,[1,M*N]),reshape(Y,[1,M*N]),reshape(z,[1,M*N]),'filled','black')
    xlabel('X')
    ylabel('Y')
    legend('Datapunten')
    title(['Benaderende 3D plot van functie van Runge tot op graad m=',num2str(m),' en n=',num2str(n),' met equidistante verdeling'])

%Niet-equidistante verdeling

    %Verdeling maken.
    x=[linspace(-1,-0.20, 5), linspace(-0.18, 0.18, 20), linspace(0.20, 1, 5)];
    y=x;
    %Puntenrooster maken.
    [X,Y]=meshgrid(x,y);
    %Functiewaarden van gegeven functie op puntenrooster.
    F=f_handle(X,Y);
    
    %Plot.
    figure()
    surf(X,Y,F)
    xlabel('X')
    ylabel('Y')
    title('Werkelijke 3D plot van functie van Runge met niet-equidistante verdeling')
    
    %Waarden opvragen.
    CoefMatrix=kkb2d(x, y, F, m, n);
    z=polyval2(CoefMatrix,x,y);%Approximatie.
    %Maximale fout bepalen.
    MaxFoutBenadering(2,1)=max(abs(F-z),[], 'all');
    %Plot.
    figure()
    surf(X,Y,z,'HandleVisibility','off')
    hold on
    scatter3(reshape(X,[1,M*N]),reshape(Y,[1,M*N]),reshape(z,[1,M*N]),'filled','black')
    xlabel('X')
    ylabel('Y')
    legend('Datapunten')
    title(['Benaderende 3D plot van functie van Runge tot op graad m=',num2str(m),' en n=',num2str(n),' met niet-equidistante verdeling'])

    %Maximale fout
disp(['Maximale fout voor equidistante verdeling: ',num2str(MaxFoutBenadering(1,1)),' en voor niet-equidistanteverdeling: ',num2str(MaxFoutBenadering(2,1)),'.'])

%Ter info de tijd die de berekening nam.  
TijdTerInfoInSeconden=toc;
disp(['Het uitvoeren van alle code duurde: ', num2str(TijdTerInfoInSeconden),' seconden.'])
%% Etna
clear
tic
%Parameters ingeven.
M=435;
N=484;
%Waarden inlezen.
F=double(imread('etna.jpg'));
%Verdeling maken.
x=linspace(-1,1,M);
y=linspace(-1,1,N);
%Puntenrooster maken.
[X,Y]=meshgrid(x,y);

%Plot.
surf(X,Y,transpose(F),'FaceLighting','gouraud','EdgeColor','none','LineStyle','none');%The 'phong' value has been removed. Use 'gouraud' instead.
colorbar
xlim([-1,1]);
ylim([-1,1]);
zlim([0,250]);
title('Etna');
figure()
imagesc(F)
colorbar
title('Bovenaanzicht Etna');

%Approximatie
m=25;
n=25;

%Waarden opvragen.
CoefMatrix=kkb2d(x, y, F, m, n);
z=polyval2(CoefMatrix,x,y);%Calculate approximation

%Plot.
figure()
surf(X,Y,z,'FaceLighting','gouraud','EdgeColor','none','LineStyle','none')%The 'phong' value has been removed. Use 'gouraud' instead.
colorbar
xlabel('X')
ylabel('Y')
title(['Benaderende 3D plot van Etna tot op graad m=',num2str(m),' en n=',num2str(n)])
figure()
imagesc(z)
colorbar
title(['Benadering bovenaanzicht Etna tot op graad m=',num2str(m),' en n=',num2str(n)])

%Ter info de tijd die de berekening nam.  
TijdTerInfoInSeconden=toc;
disp(['Het uitvoeren van alle code duurde: ', num2str(TijdTerInfoInSeconden),' seconden.'])