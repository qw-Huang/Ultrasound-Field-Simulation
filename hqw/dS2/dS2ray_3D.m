%rayleigh求三维空间中的声强和FOCUS内置rayleigh三维空间中的声场将声压转化为声强进行比较计算误差，画出来的误差图只选取了感兴趣区域
clc;
clear all;
clear all;
f0=1e6;u=1;%定义频率和法向阵速
medium = set_medium('lossless');%定义介质：单层->水，可改成多层，用函数set_layered_medium
lambda = medium.soundspeed/f0;%波长=c/f 1.5mm   lambda/6=0.25mm
k=2*pi/lambda;%波数


Roc = 5 * 2 * lambda;%ROC曲率半径
a = 5 * lambda;%注意这里a是孔径的一半
fnumber=Roc/(2*a);%所以f-number=曲率半径/孔径（2*a）
d = sqrt(Roc^2 - a^2);%理论焦点到孔径中心的距离
phi0=asin(a/Roc);
%划分网格点
xmin=-a;%观察点坐标的范围
xmax=-xmin;
ymax=xmax;
ymin=-ymax;
zDiff=0.7*d;
zmin=Roc-zDiff;
zmax=Roc+zDiff;

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
dR=lambda/6;%每个环带的弧度相等
dphi=dR/Roc;%xz平面上对应的等角度
nphi=round(phi0/dphi); %环带的个数，
phi1=dphi:dphi:phi0;
phi=phi1-dphi/2;%ri对应的xz平面的角度
ri0=Roc*sin(phi);%ri是第i环对应半径
Sm=lambda/6;%每个圆环对应的圆弧间隔频率固定
ntheta=round(2*pi*ri0/Sm);  %第i个圆环离散的弧度个数
dtheta=Sm./ri0;  %每个圆环离散的弧度
% detheta1=dtheta-dtheta/2;
theta=dthetarepet(dtheta,ntheta);%变成一个数组
ri=repelem(ri0,ntheta);
Z0=Roc-Roc*cos(phi); %每个半径为ri的圆环对应的z轴坐标
Z=repelem(Z0,ntheta);

% 求y=0时xz平面
y=0;

%rayleigh积分计算xz平面的声场  
tic
for ix=1:nx  %观察网格点x方向的坐标
    for iz=1:nz  %观察网格点z方向的坐标
        rn=sqrt(ri.^2+x(ix)^2+(Z-z(iz)).^2-2*ri.*x(ix).*cos(theta).*sign(x(ix)));%观察点到点源的距离
        dS=dR*Sm; % 求点源面积=圆环的采样间隔dR*第i个环的离散的弧度
        A=dS.*exp(-1i.*k.*rn)./rn;
        B=sum(sum(A));%对上述求得的值累加
        p=1i*medium.density*u*medium.soundspeed*k/(2*pi)*B; %乘以相关参数得到声压p
        pr(ix,iz)=squeeze(p);%得到网格点上对应的值，取模
    end
end
I_pr= acousticintensity(pr,medium.density,medium.soundspeed); %声强计算
toc

%调用focus内置rayleigh函数
xdcr = get_spherical_shell(a,Roc); %得到球面换能器
delta = [dx dy dz];%网格点步长
ps = set_coordinate_grid(delta, xmin, xmax, ymin, ymax, zmin, zmax);%网格点的划分
ndiv = 80;%积分的点数
dflag = 0;%如果=1，结果有什么不同？

tic
pre=rayleigh_cw(xdcr,ps,medium,ndiv,f0,dflag);%FOCUS自带的fnm计算
I_pre= acousticintensity(pre,medium.density,medium.soundspeed); %声强计算 
toc

%找到感兴趣的区域声强-6dB范围，画出这个范围内的error;找到最大值，返回对应的三维坐标，在最大点对应的xy平面，把感兴趣区域画出来
max_index=find_maxpoint(I_pre);%返回最大点位置坐标
x_index=max_index(1);%最大值位置x下标
y_index=max_index(2);%最大值位置y下标
z_index=max_index(3);%最大值位置z下标
I_pr_max=max(I_pr(:));%自定义rayleigh声强最大值
I_pre_max=max(I_pre(:));

%计算三维空间中两种方法的误差
I_pre_normalized=I_pre./I_pre_max;%对声强归一化
I_pr_normalized=I_pr./I_pr_max;%对声强归一化
error_I=abs(I_pr_normalized-I_pre_normalized)./I_pr_normalized;

%选择三维空间感兴趣区域误差
[index]=find(I_pr_normalized>=0.25);
error_zeros(index)=error_I(index);

%画出感兴趣区域误差
figure(1);
surf(error_zeros(:,:,z_index));
axis equal;
shading interp;
colorbar;
xlabel('x');
ylabel('y');
title('rayleigh积分3D VS FNM(FOCUS) 3D在焦点处的xy平面误差情况');

figure(2);
error_xz=squeeze(error_zeros(:,y_index,:));
surf(error_xz);
axis equal;
shading interp;
colorbar;
xlabel('z');
ylabel('x');
title('rayleigh积分3D VS FNM(FOCUS)3D在焦点处的xz平面误差情况');

figure(3);
surf(I_pre_normalized(:,:,z_index));
axis equal;
shading interp;
colorbar;
xlabel('y');
ylabel('x');
title('rayleigh(FOCUS) 3D在焦点处的xy平面');
figure(4);
surf(I_pr_normalized(:,:,z_index));
axis equal;
shading interp;
colorbar;
xlabel('y');
ylabel('x');
title('rayleigh积分3D 在焦点处的xy平面');
% 
% figure(5);
% plot(max(error_zeros(:,y_index,:)));
% xlabel('z');
% ylabel('error');
% title('rayleigh积分3D VS FNM 3D 感兴趣体间隔1mm，每个面的最大误差情况');

% 
% %轴向-6dB 网格点17-62  每隔4个点间隔1mm
% zz=17:1:61;
% for i=1:length(zz)
%     error_max(i)=max(max(error_zeros(:,:,zz(i))));
%     error_average(i)=mean(mean(error_zeros(:,:,zz(i))));
%     error_median(i)=median(error_zeros(:,:,zz(i)),'all');
% end

% figure(4);
% plot(zz,error_average);
% xlabel('z');
% ylabel('error');
% title('感兴趣体间隔1mm，每个面的平均误差情况');
% figure(5);
% plot(zz,error_median);


