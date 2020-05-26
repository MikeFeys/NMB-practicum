%% Opgave 4
clear
clc
tic

%Data inlezen
VolSymmetrischeMatrix=load('matrix.txt');
%Controleren of de matrisch hermitisch is
tf = ishermitian(VolSymmetrischeMatrix);
tfskew = ishermitian(VolSymmetrischeMatrix,'skew');
if tf
    disp('De matrix is hermitisch.')
else if tfskew
        disp('De matrix is scheef-hermitisch.')
    else
        disp('De matrix is niet hermitisch.')
    end
end
%Matrix omvormen tot hessenberg vorm.
HessenbergVorm=hess(VolSymmetrischeMatrix);

%Algemene matrix structuur weergeven van de omgezette matrix in hessenberg
%vorm.
spy(HessenbergVorm,'b',25)

TijdTerInfoInSeconden= toc;
disp(['De tijd die de berekening nam was: ',num2str(TijdTerInfoInSeconden),' seconden.'])
%% Opgave 5
clear
clc
tic

%Data inlezen.
VolSymmetrischeMatrix=load('matrix.txt');

%Eigenwaarden en residuen krijgen volgens methodes.
[EigenwaardeQRZonder,ResiduQRZonder]=qr_zonder(VolSymmetrischeMatrix);
[EigenwaardeQRWilkinson,ResiduQRWilkinson]=qr_shiftwilkinsonall(VolSymmetrischeMatrix);
[EigenwaardeQRRayleigh,ResiduQRRayleigh]=qr_shiftall(VolSymmetrischeMatrix);

%Kiezen welke eigenwaarde
eigw=1;

%Convergentie voor 1 eigenwaarde plotten met elke methode apart.
figure()
plot(1:677,log10(ResiduQRZonder(eigw,:)),':gs','MarkerSize',2)
title('QR-methode zonder shift')
xlabel('Stappen')
ylabel('Log residu')
legend('Log residuwaarden van 1 eigenwaarde')
figure()
plot(1:65,log10(ResiduQRWilkinson(eigw,:)),':ro','MarkerSize',2)
title('QR-methode met Wilkinson shift')
xlabel('Stappen')
ylabel('Log residu')
legend('Log residuwaarden van 1 eigenwaarde')
figure()
plot(1:83,log10(ResiduQRRayleigh(eigw,:)),':b*','MarkerSize',2)
title('QR-methode met Rayleigh quotiënt shift')
xlabel('Stappen')
ylabel('Log residu')
legend('Log residuwaarden van 1 eigenwaarde')
figure()

%Een aantal stappen plotten voor alle methodes
maxstappen=60;
semilogy(1:maxstappen,ResiduQRRayleigh(eigw,1:maxstappen),'b*',1:maxstappen,ResiduQRWilkinson(eigw,1:maxstappen),'-ro',1:maxstappen,ResiduQRZonder(eigw,1:maxstappen),'g--s')
%semilogy(1:30,resZonder(:,1:30),'-o', 1:size(resRay,2),resRay,'r:d', 1:size(resWilk,2),resWilk,'g--s');
title('Convergentie voor 1 eigenwaarde')
xlabel('Stappen')
ylabel('Log residu')
legend('qrRayleighShift', 'qrWilkinsonShift','qrZonder')

%Data voor fit
rangeZonder=20:677;%Eerste 19 weg
rangeRayleigh=1:81;%Laatste 2 waarden weg
rangeWilkinson=1:64;%Laatste waarde weg

dataZonder=log10(ResiduQRZonder(eigw,rangeZonder));
dataRayleigh=log10(ResiduQRRayleigh(eigw,rangeRayleigh));
dataWilkinson=log10(ResiduQRWilkinson(eigw,rangeWilkinson));

%Maak een tabel met waarden van de verschillende methodes.
format shortG
LogTabel=[ResiduQRZonder(1,:);ResiduQRRayleigh(1,:) NaN(1,677-83);ResiduQRWilkinson(1,:) NaN(1,677-65)];
TabelMetWaarden=[1,64:65,82:83, 676:677;LogTabel(:,[1,64:65,82:83, 676:677])]

TijdTerInfoInSeconden= toc;
disp(['De tijd die de berekening nam was: ',num2str(TijdTerInfoInSeconden),' seconden.'])
%% Opgave 6 
clear
clc
tic
%Tolerantie kiezen.
tol=1.e-13;

%Data inlezen.
VolSymmetrischeMatrix=load('matrix.txt');

%Jacobi methode toepassen.
[Eigenvectors,EigenwaardenJacobi] = jacobi(VolSymmetrischeMatrix,tol);%Jacobi methode moet worden aangepast om residuen terug te geven.



TijdTerInfoInSeconden= toc;
disp(['De tijd die de berekening nam was: ',num2str(TijdTerInfoInSeconden),' seconden.'])