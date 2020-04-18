%% Opgave 4
clear
clc
tic

%Data inlezen
VolSymmetrischeMatrix=load('matrix.txt');

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

%Gemiddelde residu van alle rijen opgeteld gemiddeld plotten??

%QR-methode zonder shift.
figure()
gemiddeldeResiduQRZonder=sum(ResiduQRZonder,1)/length(ResiduQRZonder(:,1));
plot(1:677,log(gemiddeldeResiduQRZonder),'b*','MarkerSize',1)
title('QR-methode zonder shift')
xlabel('Stappen')
ylabel('Log residu')
legend('Log residuwaarden')
figure()
plot(1:677,gemiddeldeResiduQRZonder,'b*','MarkerSize',1)
title('QR-methode zonder shift')
xlabel('Stappen')
ylabel('Residu')
legend('Residuwaarden')

%QR-methode met Wilkinson shift.
figure()
gemiddeldeResiduQRWilkinson=sum(ResiduQRWilkinson,1)/length(ResiduQRWilkinson(:,1));
plot(1:65,log(gemiddeldeResiduQRWilkinson),'b*','MarkerSize',3)
title('QR-methode met Wilkinson shift')
xlabel('Stappen')
ylabel('Log residu')
legend('Log residuwaarden')
figure()
plot(1:65,gemiddeldeResiduQRWilkinson,'b*','MarkerSize',3)
title('QR-methode met Wilkinson shift')
xlabel('Stappen')
ylabel('Residu')
legend('Residuwaarden')

%QR-methode met Rayleigh quoti�nt shift.
figure()
gemiddeldeResiduQRRayleigh=sum(ResiduQRRayleigh,1)/length(ResiduQRRayleigh(:,1));
plot(1:83,log(gemiddeldeResiduQRRayleigh),'b*','MarkerSize',3)
title('QR-methode met Rayleigh quoti�nt shift')
xlabel('Stappen')
ylabel('Log residu')
legend('Log residuwaarden')
figure()
plot(1:83,gemiddeldeResiduQRRayleigh,'b*','MarkerSize',3)
title('QR-methode met Rayleigh quoti�nt shift')
xlabel('Stappen')
ylabel('Residu')
legend('Residuwaarden')

%Maak een tabel met logwaarden van de verschillende methodes.
LogTabel=[1:677;log10(gemiddeldeResiduQRZonder);log10(gemiddeldeResiduQRRayleigh) NaN(1,677-83);log10(gemiddeldeResiduQRWilkinson) NaN(1,677-65)]';

%Change the table to latex format for easy use. Only useful parts are
%selected from table.
latex_table = latex(vpa(sym([LogTabel(1:5,:); LogTabel(63:66,:); LogTabel(80:84,:); LogTabel(674:677,:)]),2));
disp(['Latex formaat van de nuttige delen van de tabel: ',latex_table])

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