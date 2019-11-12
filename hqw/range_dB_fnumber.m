%�ı�f-number���̶�a���ı�R������-6dB��Χ�ľ��Դ�Сclc;clear all
clc
clear all;
f0=1e6;%    fnumber=R/(2*a);%����f-number=���ʰ뾶/�׾���2*a��
medium = set_medium('lossless');%������ʣ�����->ˮ
lambda = medium.soundspeed/f0;%����=c/f 1.5mm   lambda/6=0.25mm
fnumber=0.7:0.05:1.4;
for i=1:length(fnumber)
    a = 5 *1* 1.5e-3;%ע������a�ǿ׾���һ��
    R = 2*a*fnumber(i);%ROC���ʰ뾶
    range=range_dB_ASA(R,a,f0);
    range_radius(i)=range(1);
    range_axial(i)=range(2);
    range_z(i)=range(3);
end
figure(1);
plot(fnumber,range_radius);
title('����-6dB��Χ��С');
xlabel('f-number');
ylabel('radius');

figure(2);
plot(fnumber,range_axial);
title('����-6dB��Χ��С');
xlabel('f-number');
ylabel('axial');

figure(3);
plot(fnumber,range_z);
title('�ı�f-number������ǰ�ƾ���仯���');
xlabel('f-number');
ylabel('z');
