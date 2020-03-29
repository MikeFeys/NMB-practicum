function [T] = cheb(n,T)

if nargin < 2
    T = zeros(n+1,n+1);
end    
    
if n==0
    T(1,1:2) = [1,0];
    return;
end
if n==1
    T(1,1:2) = [1,0];
    T(2,1:2) = [0,1];
    return;
end
T = cheb(n-1,T);
T(n+1,:) = 2*circshift(T(n,:),1) - T(n-1,:);
