clc;
clear all;
%�Ƚ�ͬһ�����滻����FOCUS��ʹ��FNM������rayleigh���ֺ�������ѹ�ֲ�������,prs-pref,��ȡģ���ڱ�Ե�����ԵĲ������λ�ò�𲻴�
%����=1.5mm��R=15mm f-number=1    x����ֱ��ʣ�0.25mm   z����ֱ��ʣ�0.2598mm  ��ɢ������Դ��С������Դ����ndiv=200��
%dflagʲô���壿һ����Ϊ0���ĳ�1����
lossless = set_medium('lossless');

f0 = 1e6;
lambda = lossless.soundspeed / f0;
omega = 2 * pi * f0;
k = 2 * pi / lambda;

R =  5 * 2 * lambda;
a =  5 * lambda;
d = sqrt(R^2 - a^2);
phi0 = asin(a/R);

xdcr = get_spherical_shell(a,R);
xmin = -a;
xmax = a;
ymin = -a;
ymax = a;
z0=R;
zmin=z0;
zmax=z0;

nx = 121;
ny = 121;

dx= 2 * a / (nx - 1);
dy= 2 * a / (ny - 1);
delta = [dx dy 0];

x = xmin:dx:xmax;
y= ymin:dy:ymax;

ps = set_coordinate_grid(delta, xmin, xmax, ymin, ymax, zmin, zmax);

ndiv = 100;
dflag = 0;
tic
pref=fnm_call(xdcr,ps,lossless,ndiv,f0,dflag);
toc

tic
prs=rayleigh_cw(xdcr,ps,lossless,ndiv,f0,dflag);
RStimesphere= toc;

% plot the pressures
figure(1);
pref_normalized=abs(pref)./abs(max(pref(:)));
surf(pref_normalized);
shading interp
colorbar;
title('Fast Nearfield Method Result');
xlabel('y');
ylabel('x ');
zlabel('normalized pressure');

figure(2);
prs_normalized=abs(prs)./abs(max(prs(:)));
surf(prs_normalized);
shading interp
colorbar;
title('Rayleigh Sommerfeld Result');
xlabel('y ');
ylabel('x');
zlabel('normalized pressure');

figure(3);
plot(prs_normalized(30,:));
title('x=0 Rayleigh Sommerfeld Result');
xlabel('y');
ylabel('normalized pressure');

figure(4);
plot(prs_normalized(:,30));
title('y=0 Rayleigh Sommerfeld Result');
xlabel('x');
ylabel('normalized pressure');

error=abs(squeeze(pref - prs))./abs(squeeze(pref));%���С��1%
figure(5);
surf(error);
shading interp
colorbar;
title('error of between rayleigh and FNM');
xlabel('y ');
ylabel('x');
zlabel('pressure');

figure(6);
error0=abs(pref/abs(max(max(pref)))-prs/abs(max(max(prs))))./abs(pref/abs(max(max(pref)))); %���С��1%
surf(squeeze(error0));
shading interp
