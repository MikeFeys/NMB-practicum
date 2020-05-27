p1 = -0.1718;
p2 = 0.1718;

fit = @(x) p1*x + p2;
toplot = linspace(1,42);
hold on
plot(toplot,fit(toplot),'-r','displayname',"Lineaire fit")