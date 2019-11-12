function [theta] = dthetarepet(dtheta,ntheta)
%REPET 将环带的rn每一个都重复展开
%   此处显示详细说明
    N=1;
for i=1:length(ntheta)
    theta(N:N+ntheta(i)-1)=(dtheta(i)-dtheta(i)/2):dtheta(i):(dtheta(i)*ntheta(i)-dtheta(i)/2);
    N=N+ntheta(i);
end
[output]=theta;
end

