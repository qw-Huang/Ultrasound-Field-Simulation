%使用focus内置FNM求三维空间中的声强分布,对比和focus内置fnm+ASA的感兴趣体误差  
clc;
clear all;
clear all;
f0=1e6;u=1;%定义频率和法向阵速
medium = set_medium('lossless');%定义介质：单层->水，可改成多层，用函数set_layered_medium
lambda = medium.soundspeed/f0;%波长=c/f
k=2*pi/lambda;%波数

R = 75e-3;%ROC曲率半径
a = 30e-3;%注意这里a是孔径的一半
fnumber=R/(2*a);%所以f-number=曲率半径/孔径（2*a）
d = sqrt(R^2 - a^2);%理论焦点到孔径中心的距离
%划分网格点
xmin=-1.5*a;%观察点坐标的范围
xmax=-xmin;
ymax=xmax;
ymin=-ymax;
zDiff=0.7*d;
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

%换能器离散化为点声源  注意离散化之后点声源的大小
ntheta=100;%
dtheta=2*pi/ntheta;dr=lambda/6;
theta=dtheta:dtheta:2*pi;
r0=0:dr:a-dr;%
r=dr:dr:a;
X=cos(theta)'*r;%点源的x坐标  矩阵形式 ntheta*nr
Y=sin(theta)'*r;%点源的y坐标
Z0=R-sqrt(R*R-r.*r);%点声源三维空间中的z坐标
Z=repmat(Z0,ntheta,1); %点声源三维空间中的z坐标 repmat( A , m , n )：将向量／矩阵在垂直方向复制m次，在水平方向复制n次。

z0=zmin;  %求z=z0平面的声场分布
cg_3d = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin,zmax);%ASA求三维空间声场的网格点划分

%FNM求三维空间
xdcr = get_spherical_shell(a,R); %得到球面换能器
delta = [dx dy dz];%网格点步长
ps = set_coordinate_grid(delta, xmin, xmax, ymin, ymax, zmin, zmax);%网格点的划分
ndiv = 80;%积分的点数
dflag = 0;%如果=1，结果有什么不同？
tic
pref=fnm_call(xdcr,ps,medium,ndiv,f0,dflag);%FOCUS自带的fnm计算
I_pref= acousticintensity(pref,medium.density,medium.soundspeed); %声强计算 
toc

% FNM结合ASA计算
tic
delta2 = [dx dy 1];%网格点步长
ps2 = set_coordinate_grid(delta2, xmin, xmax, ymin, ymax, zmin, zmin);%网格点的划分
pref=fnm_call(xdcr,ps2,medium,ndiv,f0,dflag);%FOCUS自带的fnm计算
pref_asa= cw_angular_spectrum(pref,cg_3d,medium,f0,1024,'Pa');%调用内置ASA函数计三维声场分布
I_pref_asa=acousticintensity(pref_asa,medium.density,medium.soundspeed); %声强计算
toc


%找到感兴趣的区域声强-6dB范围，画出这个范围内的error;找到最大值，返回对应的三维坐标，在最大点对应的xy平面，把感兴趣区域画出来
max_index=find_maxpoint(I_pref);%返回最大点位置坐标
x_index=max_index(1);%最大值位置x下标
y_index=max_index(2);%最大值位置y下标
z_index=max_index(3);%最大值位置z下标
I_pref_asa_max=max(I_pref_asa(:));%自定义rayleigh 3D声强最大值
I_pref_max=max(I_pref(:));


error_I=abs(I_pref_asa-I_pref)./I_pref;%计算三维空间中两种方法的误差
%选择三维空间感兴趣区域误差
[index]=find(I_pref>=0.25*I_pref_max);
error_zeros(index)=error_I(index);

%画出感兴趣区域误差
figure(1);
surf(error_zeros(:,:,z_index));
shading interp;
axis equal;
xlabel('x');
ylabel('y');
title('FNM 3D VS FNM+ASA 焦点处xy平面误差情况');

figure(2);
error_xz=squeeze(error_zeros(:,y_index,:));
surf(error_xz);
shading interp;
axis equal;
xlabel('z');
ylabel('x');
title('FNM 3D VS FNM+ASA 焦点处xz平面误差情况');


figure(3);
I_pref_xz=finddB(I_pref(:,y_index,:),nx,nz);
pcolor(I_pref_xz);
axis equal
shading interp; 
colorbar;
xlabel('z');
ylabel('x');
title('FNM 焦点处xz平面声强分布情况');

figure(4);
I_pref_asa_xz=finddB(I_pref_asa(:,y_index,:),nx,nz);
pcolor(I_pref_asa_xz);
axis equal
shading interp; 
colorbar;
xlabel('z');
ylabel('x');
title('FNM+ASA 焦点处xz平面声强分布情况');
figure(6);
plot(max(error_xz));
xlabel('z');
ylabel('error');
title('FNM VS FNM+ASA 感兴趣体xz截面最大误差的变化情况');

% zz=17:1:61;
% for i=1:length(zz)
%     error_max(i)=max(max(error_zeros(:,:,zz(i))));
%     error_average(i)=mean(mean(error_zeros(:,:,zz(i))));
%     %error_median(i)=median(error_zeros(:,:,zz(i)),'all');
% end
% 
% figure(3);
% plot(zz,error_max);
% xlabel('z');
% ylabel('error');
% title('FNM 3D VS FNM+ASA 感兴趣体xy截面最大误差的变化情况');
% figure(4);
% plot(zz,error_average);
% xlabel('z');
% ylabel('error');
% title('感兴趣体xy截面平均误差的变化情况');


