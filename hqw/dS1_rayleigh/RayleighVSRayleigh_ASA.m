%使用内置Rayleigh积分求三维空间中的声强分布,对比内置rayleigh积分3D和focus内置rayleigh+ASA的感兴趣体误差  
clc;
clear all;
clear all;
f0=1e6;u=1;%定义频率和法向阵速
medium = set_medium('lossless');%定义介质：单层->水，可改成多层，用函数set_layered_medium
lambda = medium.soundspeed/f0;%波长=c/f 1.5mm   lambda/6=0.25mm
k=2*pi/lambda;%波数

R = 5 * 2 * lambda;%ROC曲率半径
a = 5 * lambda;%注意这里a是孔径的一半
fnumber=R/(2*a);%所以f-number=曲率半径/孔径（2*a）
d = sqrt(R^2 - a^2);%理论焦点到孔径中心的距离

%划分网格点
xmin=-a;%观察点坐标的范围
xmax=-xmin;
ymax=xmax;
ymin=-ymax;
zDiff=d;
zmin=R-zDiff;
zmax=R+zDiff;

dx = lambda/6; %网格点的步长
dy = lambda/6; 
dz = lambda/6;

x=xmin:dx:xmax;
y=ymin:dy:ymax;
z=zmin:dz:zmax;

nx=length(x);
ny=length(y);
nz=length(z);
error_zeros=zeros(nx,ny,nz);

% rayleigh（FOCUS）结合ASA计算
tic
xdcr = get_spherical_shell(a,R); %得到球面换能器
delta = [dx dy 0];%网格点步长
ps = set_coordinate_grid(delta, xmin, xmax, ymin, ymax, zmin, zmin);%网格点的划分
ndiv = 100;%积分的点数
dflag = 0;%如果=1，结果有什么不同？
pre=rayleigh_cw(xdcr,ps,medium,ndiv,f0,dflag);%FOCUS自带的fnm计算
z0=zmin;  %求z=z0平面的声场分布
cg_3d = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin,zmax);%ASA求三维空间声场的网格点划分
pre_asa= cw_angular_spectrum(pre,cg_3d,medium,f0,1024,'Pa');%调用内置ASA函数计三维声场分布
% I_pref_asa=acousticintensity(pref_asa,medium.density,medium.soundspeed); %声强计算
toc

%rayleigh（FOCUS）
tic
pre=rayleigh_cw(xdcr,cg_3d,medium,ndiv,f0,dflag);%FOCUS自带的rayleigh计算
toc

%找到感兴趣的区域声强-6dB范围，画出这个范围内的error;找到最大值，返回对应的三维坐标，在最大点对应的xy平面，把感兴趣区域画出来
max_index=find_maxpoint(abs(pre));%返回最大点位置坐标
x_index=max_index(1);%最大值位置x下标
y_index=max_index(2);%最大值位置y下标
z_index=max_index(3);%最大值位置z下标
pre_asa_max=max(abs(pre_asa(:)));%自定义rayleigh 3D声强最大值
pre_max=max(abs(pre(:)));

error=abs((pre_asa-pre)./pre_asa);%计算三维空间中两种方法的误差

%选择三维空间感兴趣区域误差
[index]=find(abs(pre)>=0.25*pre_max);
error_zeros(index)=error(index);

%画出感兴趣区域误差
figure(1);
axis equal;
surf(error_zeros(:,:,z_index));
shading interp
colorbar;
xlabel('x');
ylabel('y');
title('Rayleigh(FPCUS)3D VS Rayleigh(FOCUS)+ASA 焦点处xy平面误差情况');

figure(2);
axis equal;
error_xz=squeeze(error_zeros(:,y_index,:));
surf(error_xz);
shading interp
colorbar;
xlabel('z');
ylabel('x');
title('Rayleigh(FPCUS)3D VS Rayleigh(FOCUS)+ASA 焦点处xz平面误差情况');

figure(3);
plot(max(error_xz));
xlabel('z');
ylabel('error');
title('rayleigh 3D VS Rayleigh(FOCUS)+ASA 感兴趣体xy截面最大误差的变化情况');
% figure(4);
% plot(zz,error_average);
% xlabel('z');
% ylabel('error');
% title('rayleigh 3D VS Rayleigh(FOCUS)+ASA 感兴趣体xy截面平均误差的变化情况');


