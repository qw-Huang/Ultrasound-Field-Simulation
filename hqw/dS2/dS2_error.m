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
phi0=asin(a/R);

%��������ɢ��Ϊ����Դ  ע����ɢ��֮�����Դ�Ĵ�С
dR=lambda/6;%ÿ�������Ļ������
dphi=dR/R;%xzƽ���϶�Ӧ�ĵȽǶ�
nphi=round(phi0/dphi); %�����ĸ�����
dR=R*phi0/nphi;%nphiȡ��֮�����µ�������
dphi=dR/R;%����֮��dphi��ֵ
phi_after=dphi:dphi:phi0;
phi_back=0:dphi:(phi0-dphi);
phi=phi_after-dphi/2;%ri��Ӧ��xzƽ��ĽǶ�
ri0=R*sin(phi);%ri�ǵ�i����Ӧ�뾶
Sm=lambda/6;%ÿ��Բ����Ӧ��Բ�����Ƶ�ʹ̶�
ntheta=round(2*pi*ri0/Sm);  %��i��Բ����ɢ�Ļ��ȸ���
dtheta=2*pi./ntheta; %����ȡ�����������µ���ÿһ����dtheta
Sm0=dtheta.*ri0; %����ÿһ����ɢ�Ļ���
Sm=repelem(Sm0,ntheta);
theta=dthetarepet(dtheta,ntheta);%���һ������
ri=repelem(ri0,ntheta);
Z0=R-R*cos(phi); %ÿ���뾶Ϊri��Բ����Ӧ��z������
Z=repelem(Z0,ntheta);

%rayleigh���ּ���xzƽ�������  
tic
dS=dR.*Sm; % ���Դ���=Բ���Ĳ������dR*��i��������ɢ�Ļ���
toc

%�Ƚ�dS���ۼ�
S_discrete=sum(dS);%��ɢ֮���ȫ������ۼ�
S_theory=2*pi*R*(R-d);%���ۼ���������������������
dSi=Sm0.*dR.*ntheta;%��ÿһ��������ۼӣ�Ҳ���ǰ�ÿһ�У�70�������ۼ�
ri_back=R*sin(phi_back);%ri�ǵ�i���̻��˶�Ӧ�뾶
ri_after=R*sin(phi_after);%ri�ǵ�i����=�����˶�Ӧ�뾶
di_back=sqrt(R^2-ri_back.^2);%��0��ʼ��Ϊ��i�㻷��d
di_after=sqrt(R^2-ri_after.^2);%��dr��ʼ��Ϊ��i�㻷��d
dSi_theory=2*pi*R*(R-di_after)-2*pi*R*(R-di_back);%��i���Ļ������
error_dSi=abs(dSi-dSi_theory)./dSi;

%��ͼ
figure(4);
histogram(dS,10);
title('��Դ���dS�ֲ�(R=150mm ,a=75mm)');

figure(5);
plot(error_dSi);
xlabel('ith ring ');
ylabel('error ');
title('the error of ith ring (R=150mm, a=75mm)');
