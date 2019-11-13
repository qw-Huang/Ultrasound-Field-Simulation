%������FNM��������xzƽ����ѹ�ֲ����õ�ƽ������;�����ѹ-3dB������ǿ-6dB����С
clc;
clear all;
lossless = set_medium('lossless');
f = 1e6;
lambda = lossless.soundspeed / f;
omega = 2 * pi * f;
k = 2 * pi / lambda;

R = 6 * 1.1*lambda;
a = 6* lambda;
d = sqrt(R^2 - a^2);
phi0 = asin(a/R);
dBr=1.1*lambda;%��֤����-6dB��Χ=f-number*�������Ƿ�����ǿ��-6dB��

xdcr = get_spherical_shell(a,R);% ���λ�����
%���������
xmin = -a;
xmax = a;
ymin = 0;
ymax = 0;
zmin = R-d;
zmax = 1.5*d + R;

nz = 201; % ok to sample the origin
nx = 101;

if nx > 1,
    dx = 2 * a / (nx - 1);
else
    dx = 0;
end

if nz > 1,
    dz = 2.5*d/ (nz - 1);
else
    dz = 0;
end
delta = [dx 0 dz];

x = xmin:dx:xmax;
y = ymin:dx:ymax;
z = zmin:dz:zmax;

ps = set_coordinate_grid(delta, xmin, xmax, ymin, ymax, zmin, zmax);

ndiv = 200;
dflag = 0;
tic
pref=fnm_call(xdcr,ps,lossless,ndiv,f,dflag);
toc

figure(1)
mesh(z*1000, x*1000, abs(squeeze(pref)))
xlabel('z (mm)');
ylabel('x (mm)');
zlabel('normalized pressure');

p= abs(squeeze(pref));
[Pm1,Ix]=max(p);%�ҵ���ѹÿһ�����ֵ�����ض�Ӧ����
[Pm,Iy]=max(Pm1());%�ҵ������ֵ��ѹ��һ�У���������
Ij=Iy;
Ii=Ix(1,Ij);
Pmax=p(Ii,Ij);%�ҵ���ѹ���ֵPmax

pa=p(Ii,:);%xΪ0����һ�У�Ҳ����������ѹ�ֲ�
figure(2);
plot(z*1000,pa);
xlabel('z (mm)');
ylabel('axial pressure');

co=find(pa>0.5*Pmax);%����ѹ������У�ѡ����ѹ>0.5*Pmax����
co1=min(co);co2=max(co);%�ҵ���ѹ>0.5*Pmax����С���к���������
range_c=(co2-co1)*dz;%�и����⣬�õ��ķ�Χ�����˽�������0.25��pmax����ֵ�����ǲ�������Ҫ�ģ�ϣ���õ�������-6dB�ķ�Χ
pr=p(:,Ij);%xΪ0����һ�У�Ҳ����������ѹ�ֲ�

figure(3);
plot(x*1000,pr);
xlabel('x (mm)');
ylabel('radial pressure ');

ro=find(pr>0.5*Pmax);%����ѹ������У�ѡ����ѹ>0.5*Pmax����
ro1=min(ro);ro2=max(ro);%�ҵ���ѹ>0.5*Pmax����С���к���������
range_r=(ro2-ro1)*dx;
