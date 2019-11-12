clc;
%�Ƚ�ͬһ�����滻����FOCUS��ʹ��FNM������rayleigh���ֺ�������ѹ�ֲ�������,prs-pref,��ȡģ���ڱ�Ե�����ԵĲ������λ�ò�𲻴�
%����=1.5mm��R=15mm f-number=1    x����ֱ��ʣ�0.25mm   z����ֱ��ʣ�0.2598mm  ��ɢ������Դ��С������Դ����ndiv=200��
%RStimesphere�õ�����1:200���������ѹ��ʱ�䣬���ǲ���ֱ�ӳ��Ե�����Ϊÿ����ļ���ʱ��
%dflagʲô���壿һ����Ϊ0���ĳ�1����
clear all;
lossless = set_medium('lossless');

f0 = 1e6;
lambda = lossless.soundspeed / f0;
omega = 2 * pi * f0;
k = 2 * pi / lambda;

R = 5 * 2 * lambda;
a = 5 * lambda;
d = sqrt(R^2 - a^2);
phi0 = asin(a/R);

xdcr = get_spherical_shell(a,R);
xmin = -a;
xmax = a;
y0=-0*a;
ymin = y0;
ymax = y0;
zmin = -0.7*d + R;
zmax = 0.7*d + R;

nz = 101; % ok to sample the origin
nx = 61;
dx=(2*a)/(nx-1);
dz=(2*d)/(nz-1);
delta = [dx 0 dz];

x = xmin:dx:xmax;
z = zmin:dz:zmax;

ps = set_coordinate_grid(delta, xmin, xmax, ymin, ymax, zmin, zmax);

ndiv = 100;
dflag = 0;
tic
pref=fnm_call(xdcr,ps,lossless,ndiv,f0,dflag);
toc

RSerrorsphere = zeros(size(ndiv));

    tic
    prs=rayleigh_cw(xdcr,ps,lossless,ndiv,f0,dflag);
    toc




% plot the pressures
if nx > 1 & nz > 1,
    figure(1);
    p1=abs(squeeze(pref))/max(max(abs(squeeze(pref))));
    surf(p1);
    title('Fast Nearfield Method Result');
    xlabel('z (mm)');
    ylabel('x (mm)');
    zlabel('normalized pressure');

    figure(2);
    p2=abs(squeeze(prs))/max(max(abs(squeeze(prs))));
    surf(p2);
    title('Rayleigh Sommerfeld Result');
    xlabel('z (mm)');
    ylabel('x (mm)');
    zlabel('normalized pressure');
end

error=abs(squeeze(pref - prs))./abs(squeeze(pref));%���С��1%
figure(3);
surf(error);
title('error of between rayleigh and FNM');
xlabel('z (mm)');
ylabel('x (mm)');
zlabel('pressure');

figure(4);
error0=abs(pref/abs(max(max(pref)))-prs/abs(max(max(prs))))./abs(pref/abs(max(max(pref)))); %���С��1%
mesh(z*1000,x*1000,squeeze(error0));
datestr(now)