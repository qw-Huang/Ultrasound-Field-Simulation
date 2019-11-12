clc;
clear all;

f0=1e6;u=1;%����Ƶ�ʺͷ�������
medium = set_medium('lossless');%������ʣ����㣺ˮ�������Ըĳɶ��set_layered_medium
lambda = medium.soundspeed/f0;%����=c/f
k=2*pi/lambda;%����

R = 6 * 5 * 2 * lambda;%ROC���ʰ뾶
a = 6 * 5 * lambda;%ע�������a�ǿ׾���һ��
fnumber=R/(2*a);%����f-number=���ʰ뾶/�׾���2*a��
d = sqrt(R^2 - a^2);%���۽��㵽�׾����ĵľ���
phi0=asin(a/R);

%���������
xmin=-a;%�۲������ķ�Χ
xmax=-xmin;
ymax=0;
ymin=0;
zDiff=0.7*d;
zmin=R-zDiff;
zmax=R+zDiff;

dx = lambda/6; %�����Ĳ��� 
dz = lambda/6;

x=xmin:dx:xmax;
z=zmin:dz:zmax;

nx=length(x);
nz=length(z);
error_zeros_xz=zeros(nx,nz);
% ��y=0ʱxzƽ��
y=0;

%����ת��
r=sqrt(x.^2+y.^2);
o_theta=acos(x./r);


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
for ix=1:nx  %�۲������x���������
    for iz=1:nz  %�۲������z���������
       if x(ix)==0
           rn=sqrt(ri.^2+r(ix).^2+(Z-z(iz)).^2-2*ri.*r(ix));%�۲�㵽��Դ�ľ���
       else
           rn=sqrt(ri.^2+r(ix).^2+(Z-z(iz)).^2-2*ri.*r(ix).*cos(theta-o_theta(ix)));%�۲�㵽��Դ�ľ���
       end
        dS=dR.*Sm; % ���Դ���=Բ���Ĳ������dR*��i��������ɢ�Ļ���
        A=dS.*exp(-1i.*k.*rn)./rn;
        B=sum(sum(A));%��������õ�ֵ�ۼ�
        p=1i*medium.density*u*medium.soundspeed*k/(2*pi)*B; %������ز����õ���ѹp
        pr(ix,iz)=squeeze(p);%�õ�������϶�Ӧ��ֵ��ȡģ
    end
end
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

%����focus����rayleigh����rayleigh_cw()
omega = 2 * pi * f0;
phi0 = asin(a/R);
xdcr = get_spherical_shell(a,R); %�õ����滻����

delta = [dx 0 dz];%����㲽��
ps = set_coordinate_grid(delta, xmin, xmax, ymin, ymax, zmin, zmax);%�����Ļ���
ndiv = 100;%���ֵĵ���
dflag = 0;%���=1�������ʲô��ͬ��

tic
prs=rayleigh_cw(xdcr,ps,medium,ndiv,f0);%FOCUS�Դ���rayleigh����
toc

I_prs=acousticintensity(prs,medium.density,medium.soundspeed); 
I_prs_nor=squeeze(I_prs./max(I_prs(:)));
I_pr=acousticintensity(pr,medium.density,medium.soundspeed); 
I_pr_nor=squeeze(I_pr./max(I_pr(:)));
error=abs(I_pr_nor-I_prs_nor)./I_pr_nor;


%��ͼ
figure(1);
surf(z*1000,x*1000,I_pr_nor); %������Ϊz��������Ϊx����ɫ����p1�Ĵ�С  ��Ϊʲô�����������Ķ�Ӧ��ϵ��
shading interp
title('Rayleigh ');
axis equal;%��������ı�����ͬ
colorbar
xlabel('z ');
ylabel('x ');
zlabel('normalized pressure');

figure(2);
surf(z*1000, x*1000, I_prs_nor);%��prs��һ��   Ϊʲô��Ӧ��������z��ǰ�棿
shading interp;
colorbar;
axis equal;
title('rayleigh(FOCUS)');
xlabel('z ');
ylabel('x ');
zlabel('normalized pressure');

%ѡ�����Ȥ�������
[index]=find(I_prs_nor>=0.25);
error_zeros_xz(index)=error(index);

%�Ա����ַ����������Ľ��
figure(3);
surf(z*1000,x*1000,error_zeros_xz); 
shading interp %
axis equal;
colorbar  %����ɫ��
xlabel('z ');
ylabel('x');
title('error(R=120mm,a=60mm)');

figure(4);
% plot(z*1000,max(error_zeros_xz));
% datestr(now)
 histogram(dS,10);


figure(5);
plot(error_dSi);