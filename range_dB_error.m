%����range_dB function �õ���Ӧƽ�������;���-6dB�Ĵ�С��������Ϊ���ı�f-number���õ���ͬ��-6dB
%range��ֵ��������ֵ���Ƚϣ����������
clc;clear all;
f0=1e6;%    fnumber=R/(2*a);%����f-number=���ʰ뾶/�׾���2*a��
medium = set_medium('lossless');%������ʣ�����->ˮ
lambda = medium.soundspeed/f0;%����=c/f 1.5mm   lambda/6=0.25mm
fnumber=0.7:0.05:1.35;
for i=1:length(fnumber)
    a = 5 * 1.5e-3;%ע������a�ǿ׾���һ��
    R = 2*a*fnumber(i);%ROC���ʰ뾶
    range=range_dB_ASA(R,a,f0);
    range_radius(i)=range(1);
    range_axial(i)=range(2);
end
radius_ref=lambda.*fnumber.^2;%����-6dB�ο�ֵ
axial_ref=7*lambda.*fnumber.^2;%����-6dB
error_radius=abs(range_radius-radius_ref)./radius_ref;
error_axial=abs(range_axial-axial_ref)./axial_ref;

figure(1);
plot(fnumber,error_radius);
title('�ı�f-number������-6dB��Χ��С');
xlabel('f-number');
ylabel('radius error');

figure(2);
plot(fnumber,error_axial);
title('�ı�f-number������-6dB��Χ��С');
xlabel('f-number');
ylabel('axial error');

