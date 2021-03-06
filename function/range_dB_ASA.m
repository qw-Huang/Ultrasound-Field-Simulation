function [output]=range_dB_ASA(R,a,f0)
%rayleigh+ASA计算3D声强分布，求焦点径向、轴向-6dB大小的值和焦点前移的距离
%rayleigh+ASA
%clc;clear;close all;
f0=1e6;
u=1;%定义频率和法向阵速
medium = set_medium('lossless');%定义介质：单层->水
lambda = medium.soundspeed/f0;%波长=c/f 1.5mm   lambda/6=0.25mm
k=2*pi/lambda;%波数

%fnumber=1.15;
%a = 5 * 1.5e-3;%注意这里a是孔径的一半
%R = 2*a*fnumber;%ROC曲率半径
d = sqrt(R^2 - a^2);%理论焦点到孔径中心的距离

%划分网格点
xmin=-a;%观察点坐标的范围
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

%换能器离散化为点声源  注意离散化之后点声源的大小
dr=lambda/6;%球面分割成很多个圆环，对应的半径是r，dr是半径增加的步长
r_back=0:dr:a-dr;%点声源前一段弧长对应的r
r_after=dr:dr:a;%点声源后一段弧长对应的r
r=r_after-dr/2;%第i个环带对应的中心点的r
Sm=lambda/6;%最外层环带（2pi*a）离散化后对应的dS的弧长小于lambda/6
ntheta=round(2*pi*(a-dr/2)/Sm);%最外层一个环带的划分点数,取整
dtheta=2*pi./ntheta;%根据取整重新调整每个点声源的对应弧度
theta_after=dtheta:dtheta:2*pi;%每个环带离散成多个点对应的弧度数组
theta_back=0:dtheta:(2*pi-dtheta);%每个环带离散成多个点对应的弧度数组
theta=theta_after-dtheta/2;%第i个环带离散成点声源，dS中间的点对应的弧度数组
X=cos(theta)'*r;%点源的x坐标  矩阵形式 ntheta*nr
Y=sin(theta)'*r;%点源的y坐标
Z_ring=R-sqrt(R*R-r.*r);%点声源三维空间中的z坐标
Z=repmat(Z_ring,ntheta,1); % repmat( A , m , n )：将向量／矩阵在垂直方向复制m次，在水平方向复制n次。

% z0=zmin;  %求z=z0平面的声场分布
cg_3d = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin,zmax);%ASA求三维空间声场的网格点划分

%自定义rayleigh积分计算xy面的声场  
tic
for ix=1:nx  %观察网格点x方向的坐标划分
    for iy=1:ny %观察网格点y方向的坐标划分
            rn=sqrt((X-x(ix)).^2+(Y-y(iy)).^2+(Z-zmin).^2);%观察点到点源的距离
            dS_ring=r.*dtheta.*R.*(asin(r_after./R)-asin(r_back./R));%每一环带离散的dS的大小，第i环离散的点声源面积dS=第i环的离散弧长*dr对应的短弧
            dS=repmat(dS_ring,ntheta,1); % repmat( A , m , n )：将向量／矩阵在垂直方向复制m次，在水平方向复制n次。
            A=dS.*exp(-1i.*k.*rn)./rn;
            B=sum(sum(A));%对上述求得的值累加
            pr(ix,iy)=1i*medium.density*u*medium.soundspeed*k/(2*pi)*B; %乘以相关参数得到声压p        
    end
end

%计算3D声场分布
p_asa = cw_angular_spectrum(pr,cg_3d,medium,f0,1024,'Pa');%调用内置ASA函数计三维声场分布
I_asa=acousticintensity(p_asa,medium.density,medium.soundspeed); %声强计算
toc

%找到最大点对应的索引坐标
max_index=find_maxpoint(I_asa);%返回最大点位置坐标
x_index=max_index(1);%最大值位置x下标
y_index=max_index(2);%最大值位置y下标
z_index=max_index(3);%最大值位置z下标

%计算前移距离
z_center=fix((nz-1)/2+1);%R对应的z的下标
range_value_z=(z_center-z_index)*dz;%前移距离

%，找到最大点轴向径向对应的切面，把符合-6dB范围的声强放到一个新的数组里
I_asa_xy=squeeze(I_asa(:,:,z_index));%声强最大点对应的xy面
I_asa_xz=squeeze(I_asa(:,y_index,:));%声强最大点对应的xz面
I_range_xy=finddB(I_asa_xy,nx,ny);%xy平面找到径向-6dB范围
I_range_xz=finddB(I_asa_xz,nx,nz);%xz平面找到径向-6dB范围

%xy返回有数值的最远的两个点的列数的差，相减乘以步长就是-6dB范围的距离
range_point_xy=find_range_point(I_range_xy,nx,ny);  %可能存在数据越界，如果在边界处数据不为0，函数就会出现问题
range_point_xz=find_range_point(I_range_xz,nx,nz);

range_value_xy=range_point_xy*dy;%数值=点数*步长
range_value_xz=range_point_xz*dz;

output=[range_value_xy,range_value_xz,range_value_z];
end