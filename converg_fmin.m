function [max_fout] = converg_fminbnd(f_handle, max_ord)

for n = 1:max_ord  
a = chebcoeff(f_handle,n);
T = cheb(n);
c = poly(a,T);
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
end
figure()
plot(1:n,max_fout)
end