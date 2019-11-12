%function [output]=rayleigh_ASA(R,a,u,f0)
%输入R,a,u,f0,输出中间开孔的3D声强分布，使用自定义rayleigh积分+ASA计算声强
clc;clear all;
f0=1e6;u=1;%定义频率和法向阵速

a = 5 * 1.5e-3;%注意这里a是孔径的一半
R = 2*a;%ROC曲率半径
rhole=0.004;
medium = set_medium('lossless');%定义介质：单层->水
lambda = medium.soundspeed/f0;%波长=c/f 1.5mm   lambda/6=0.25mm
k=2*pi/lambda;%波数
%fnumber=1;
d = sqrt(R^2 - a^2);%理论焦点到孔径中心的距离
%划分网格点
xmin=-0.7*a;%观察点坐标的范围
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
ntheta=70;%
dtheta=2*pi/ntheta;dr=lambda/6;
theta=dtheta:dtheta:2*pi;
r0=rhole:dr:a-dr;%
r=(rhole+dr):dr:a;
X=cos(theta)'*r;%点源的x坐标  矩阵形式 ntheta*nr
Y=sin(theta)'*r;%点源的y坐标
Z0=R-sqrt(R*R-r.*r);%点声源三维空间中的z坐标
Z=repmat(Z0,ntheta,1); %点声源三维空间中的z坐标 repmat( A , m , n )：将向量／矩阵在垂直方向复制m次，在水平方向复制n次。

z0=zmin;  %求z=z0平面的声场分布
cg_3d = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin,zmax);%ASA求三维空间声场的网格点划分

%自定义rayleigh积分计算xy面的声场  
tic
for ix=1:nx  %观察网格点x方向的坐标划分
    for iy=1:ny %观察网格点y方向的坐标划分
           rn=sqrt((X-x(ix)).^2+(Y-y(iy)).^2+(Z-z0).^2);%观察点到点源的距离
           dS0=r.*dtheta.*R.*(asin(r./R)-asin(r0./R));%点源的面积近似为矩形,面积=长（圆弧）×宽（圆弧）
           dS=repmat(dS0,ntheta,1); % repmat( A , m , n )：将向量／矩阵在垂直方向复制m次，在水平方向复制n次。
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
I_xz=squeeze(I_asa(:,y_index,:));
mesh(I_xz);
%end