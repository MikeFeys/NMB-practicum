
p1 =  -0.01646;
p2 =   0.08211;
p3 =   0.1883;
fit = @(x) p1*x.^2 + p2*x + p3;
toplot = linspace(1,30)
hold on
plot(toplot,fit(toplot),'-r','displayname',"Kwadratische fit")