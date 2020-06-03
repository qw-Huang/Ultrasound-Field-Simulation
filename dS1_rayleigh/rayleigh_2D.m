%使用dS1点声源划分方法，在2d平面xz面，归一化对比自定义rayleigh积分和focus自带
%的rayleigh积分，对比声强分布的误差
% clc;
clear all;
f0=0.8e6;%定义频率和法向阵速
P=100;
medium = set_medium('muscle');%定义介质（单层：水），可以改成多层set_layered_medium
lambda = medium.soundspeed/f0;%波长=c/f
dBperNeper = 20 * log10(exp(1));
% dBperNeper =8.6886;
attenuationNeperspermeter=medium.attenuationdBcmMHz/dBperNeper*100*(f0/1e6)^1.1;
kk=2*pi/lambda-j*attenuationNeperspermeter;%波数
k=2*pi/lambda;%波数

% k=2*pi/lambda-j*medium.attenuationdBcmMHz;%波数

R = 150e-3;%ROC曲率半径
a=60e-3;
fnumber=R/(2*a);%所以f-number=曲率半径/孔径（2*a）
d = sqrt(R^2 - a^2);%理论焦点到孔径中心的距离
u=normal_velocity(P,R,a,0,medium.density,medium.soundspeed);
% u=1;

%划分网格点
xmin=-1.5*a;%观察点坐标的范围
xmax=-xmin;
ymax=0;
ymin=0;

zmin=140e-3;
% zmin=R-zDiff;
zmax=150e-3;

dx = lambda/6; %网格点的步长
dz = lambda/6;

x=xmin:dx:xmax;%网格点的分布
z=zmin:dz:zmax;

nx=length(x);%网格点的点数
nz=length(z);

error_zeros_xz=zeros(nx,nz); %定义全零误差矩阵

%换能器离散化为点声源  注意离散化之后点声源的大小
dr=lambda/6;%球面分割成很多个圆环，对应的半径是r，dr是半径增加的步长
r_back=0:dr:a-dr;%点声源前一段弧长对应的r
r_after=dr:dr:a;%点声源后一段弧长对应的r
r=r_after-dr/2;%第i个环带对应的中心点的r
Sm=lambda/6;%中间环带离散化后对应的dS的弧长等于lambda/6
median=round(length(r)/2);%取中间环带的索引
ntheta=round(2*pi*r(median)/Sm);%中间层一个环带的划分点数,取整
dtheta=2*pi./ntheta;%根据取整重新调整每个点声源的对应弧度
theta_after=dtheta:dtheta:2*pi;%每个环带离散成多个点对应的弧度数组
theta_back=0:dtheta:(2*pi-dtheta);%每个环带离散成多个点对应的弧度数组
theta=theta_after-dtheta/2;%第i个环带离散成点声源，dS中间的点对应的弧度数组
X=cos(theta)'*r;%点源的x坐标  矩阵形式 ntheta*nr
Y=sin(theta)'*r;%点源的y坐标
Z0=R-sqrt(R*R-r.*r);%点声源三维空间中的z坐标
Z=repmat(Z0,ntheta,1); % repmat( A , m , n )：将向量／矩阵在垂直方向复制m次，在水平方向复制n次。

% 求y=0时xz平面
y=0;

%rayleigh积分计算xz平面的声场  
tic
for ix=1:nx  %观察网格点x方向的坐标
    for iz=1:nz  %观察网格点z方向的坐标
        rn=sqrt((X-x(ix)).^2+(Y-y).^2+(Z-z(iz)).^2);%观察点到点源的距离
        dS_ring=r.*dtheta.*R.*(asin(r_after./R)-asin(r_back./R));%每一环带离散的dS的大小，第i环离散的点声源面积dS=第i环的离散弧长*dr对应的短弧
        dS=repmat(dS_ring,ntheta,1); % repmat( A , m , n )：将向量／矩阵在垂直方向复制m次，在水平方向复制n次。
        A=dS.*exp(-1i.*kk.*rn)./rn;
        B=sum(sum(A));%对上述求得的值累加 
        pr(ix,iz)=1i*medium.density*u*medium.soundspeed*k/(2*pi)*B; %乘以相关参数得到声压p
    end
end
toc

pr_max=max(abs(pr(:)));
pr_abs=abs(pr);
focus_index_1=find(pr_abs==pr_max);
s=size(pr_abs);
[x_index,z_index]=ind2sub(s,focus_index_1);%将最大值单下标转为三维多下标
focus_index=z_index;
% 
% %比较dS的累加和解析面积
% S_discrete=sum(sum(dS));%离散之后的全部面积累加
% S_theory=2*pi*R*(R-d);%理论计算球面面积，即解析面积
% dSi=sum(dS);%求每一环面积的累加，也就是把每一列（70个数）累加
% di_back=sqrt(R^2-r_back.^2);%从0开始作为第i层环的d
% di_after=sqrt(R^2-r_after.^2);%从dr开始作为第i层环的d
% dSi_theory=2*pi*R*(R-di_after)-2*pi*R*(R-di_back);%第i环的环带面积
% error_dSi=abs(dSi-dSi_theory)./dSi;%求解析解环带面积和数值解环带面积的误差
% 
% %调用focus内置rayleigh函数rayleigh_cw()
% omega = 2 * pi * f0;
% phi0 = asin(a/R);
% xdcr = get_spherical_shell(a,R); %得到球面换能器
% delta = [dx 0 dz];%网格点步长
% ps = set_coordinate_grid(delta, xmin, xmax, ymin, ymax, zmin, zmax);%网格点的划分
% ndiv = 100;%积分的点数
% dflag = 0;%如果=1，结果有什么不同？
% 
% tic
% % prs=rayleigh_cw(xdcr,ps,medium,ndiv,f0);%FOCUS自带的rayleigh计算
% prs=cw_pressure(xdcr,ps,medium,ndiv,f0);%FOCUS自带的fnm计算
% 
% toc
% 
% %声压转化为声场计算
% I_prs=acousticintensity(prs,medium.density,medium.soundspeed);
% I_prs_nor=squeeze(I_prs./max(I_prs(:)));
% I_pr=acousticintensity(pr,medium.density,medium.soundspeed); 
% I_pr_nor=I_pr./max(I_pr(:));
% error=abs(I_pr_nor-I_prs_nor)./I_pr_nor;
% 
% %选择感兴趣区域误差
% [index]=find(I_prs_nor>=0.25);
% error_zeros_xz(index)=error(index);

% figure(3)
% surf(z*1000,x*1000,abs(pr));
% shading interp;
% colorbar;
% axis equal;
