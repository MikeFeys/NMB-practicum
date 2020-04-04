%input:
%f is a function handle
% n is the maximal order of approach
function [max_fout] = converg_fminbnd(f_handle, max_ord)

max_fout = zeros(1,max_ord);
for n = 1:max_ord
    
a = chebcoeff(f_handle,n);
T = cheb(n);
c = poly(a,T);
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

end