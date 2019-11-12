%ʹ��dS1����Դ���ַ�������2dƽ��xz�棬��һ���Ա��Զ���rayleigh���ֺ�focus�Դ�
%��rayleigh���֣��Ա���ǿ�ֲ������
clc;
clear all;
f0=1e6;u=1;%����Ƶ�ʺͷ�������
medium = set_medium('lossless');%������ʣ����㣺ˮ�������Ըĳɶ��set_layered_medium
lambda = medium.soundspeed/f0;%����=c/f
k=2*pi/lambda;%����

R = 10 * 5 * 2 * lambda;%ROC���ʰ뾶
a = 10 * 5 * lambda;%ע�������a�ǿ׾���һ��
fnumber=R/(2*a);%����f-number=���ʰ뾶/�׾���2*a��
d = sqrt(R^2 - a^2);%���۽��㵽�׾����ĵľ���

%��������ɢ��Ϊ����Դ  ע����ɢ��֮�����Դ�Ĵ�С
dr=lambda/6;%����ָ�ɺܶ��Բ������Ӧ�İ뾶��r��dr�ǰ뾶���ӵĲ���
r_back=0:dr:a-dr;%����Դǰһ�λ�����Ӧ��r
r_after=dr:dr:a;%����Դ��һ�λ�����Ӧ��r
r=r_after-dr/2;%��i��������Ӧ�����ĵ��r
Sm=lambda/6;%����㻷����2pi*a����ɢ�����Ӧ��dS�Ļ���С��lambda/6
median=length(r)/2+1;
ntheta=round(2*pi*r(median)/Sm);%�����һ�������Ļ��ֵ���,ȡ��
dtheta=2*pi./ntheta;%����ȡ�����µ���ÿ������Դ�Ķ�Ӧ����
theta_after=dtheta:dtheta:2*pi;%ÿ��������ɢ�ɶ�����Ӧ�Ļ�������
theta_back=0:dtheta:(2*pi-dtheta);%ÿ��������ɢ�ɶ�����Ӧ�Ļ�������
theta=theta_after-dtheta/2;%��i��������ɢ�ɵ���Դ��dS�м�ĵ��Ӧ�Ļ�������
X=cos(theta)'*r;%��Դ��x����  ������ʽ ntheta*nr
Y=sin(theta)'*r;%��Դ��y����
Z0=R-sqrt(R*R-r.*r);%����Դ��ά�ռ��е�z����
Z=repmat(Z0,ntheta,1); % repmat( A , m , n )���������������ڴ�ֱ������m�Σ���ˮƽ������n�Ρ�

%���㻮��֮���dS
tic
dS_ring=r.*dtheta.*R.*(asin(r_after./R)-asin(r_back./R));%ÿһ������ɢ��dS�Ĵ�С����i����ɢ�ĵ���Դ���dS=��i������ɢ����*dr��Ӧ�Ķ̻�
dS=repmat(dS_ring,ntheta,1); % repmat( A , m , n )���������������ڴ�ֱ������m�Σ���ˮƽ������n�Ρ�
toc

%�Ƚ�dS���ۼӺͽ������
S_discrete=sum(sum(dS));%��ɢ֮���ȫ������ۼ�
S_theory=2*pi*R*(R-d);%���ۼ���������������������
dSi=sum(dS);%��ÿһ��������ۼӣ�Ҳ���ǰ�ÿһ�У�70�������ۼ�
di_back=sqrt(R^2-r_back.^2);%��0��ʼ��Ϊ��i�㻷��d
di_after=sqrt(R^2-r_after.^2);%��dr��ʼ��Ϊ��i�㻷��d
dSi_theory=2*pi*R*(R-di_after)-2*pi*R*(R-di_back);%��i���Ļ������
error_dSi=abs(dSi-dSi_theory)./dSi;%������⻷���������ֵ�⻷����������

%��ͼ
figure(4);
histogram(dS_ring,10);
title('��Դ���dS�ֲ�(R=150mm ,a=75mm)');

figure(5);
plot(error_dSi);
xlabel('ith ring ');
ylabel('error ');
title('the error of ith ring (R=150mm, a=75mm)');
