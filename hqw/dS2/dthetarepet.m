function [theta] = dthetarepet(dtheta,ntheta)
%REPET ��������rnÿһ�����ظ�չ��
%   �˴���ʾ��ϸ˵��
    N=1;
for i=1:length(ntheta)
    theta(N:N+ntheta(i)-1)=(dtheta(i)-dtheta(i)/2):dtheta(i):(dtheta(i)*ntheta(i)-dtheta(i)/2);
    N=N+ntheta(i);
end
[output]=theta;
end

