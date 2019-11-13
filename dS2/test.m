clear all;
dtheta=[2 3 4];

ntheta =[3 4 5];
    N=1;
for i=1:length(ntheta)
    theta(N:N+ntheta(i)-1)=dtheta(i):dtheta(i):dtheta(i)*ntheta(i);
    N=N+ntheta(i);
end