function [D]=cheb(N,a,b)
% Input:
% N - discretisation index (take N+1 Chebyshev nodes)
% a - left endpoint of  the interval
% b - right endpoint of the interval
% Output:
% D - (N+1) x (N+1) differentiation matrix
% see Trefethen

if N==0
    x=1;
    D=0;
    return
end
p=pi*(0:N)'/N;
x=((b-a)*cos(p)+b+a)/2;

c=[2;ones(N-1,1);2].*(-1).^(0:N)';
X=x(:,ones(1,N+1)); %X=repmat(x,1,N+1);
dX=X-X';
D=(c*(1./c)')./(dX+(eye(N+1)));
D=D-diag(sum(D,2)); %D=D-diag(sum(D'));