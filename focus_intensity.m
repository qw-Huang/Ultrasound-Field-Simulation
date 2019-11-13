
%FOCUS_INTENSITY ���뷨�����٣�������㴦����ǿ
%   �˴���ʾ��ϸ˵��
clc
clear all;
f0=1e6;%    fnumber=R/(2*a);%����f-number=���ʰ뾶/�׾���2*a��
medium = set_medium('lossless');%������ʣ�����->ˮ
lambda = medium.soundspeed/f0;%����=c/f 1.5mm   lambda/6=0.25mm
fnumber=0.7:0.05:1.4;
P=10;
for i=1:length(fnumber)
    a = 5 *10* 1.5e-3;%ע������a�ǿ׾���һ��
    R = 2*a*fnumber(i);%ROC���ʰ뾶
    u=normal_velocity(P,R,a,medium.density,medium.soundspeed);
    I=rayleigh_ASA(R,a,u,f0);
    I_focus(i)=max(I(:));
end
figure(1);
plot(fnumber,I_focus);
title('�ı�f-number��������ǿ�ı仯');
xlabel('f-number');
ylabel('focus intensity ��Pa��');


