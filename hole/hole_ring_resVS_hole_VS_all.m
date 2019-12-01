%使用dS1点声源划分方法，在2d平面xz面，加开环形孔的换能器声压计算，比较换能器完整部分、开孔部分、剩余部分的声压幅值和相位
clc;
close all;
clear all;
f0=1e6;%定义频率和声功率
P=100;
medium = set_medium('lossless');%定义介质（单层：水），可以改成多层set_layered_medium
lambda = medium.soundspeed/f0;%波长=c/f
k=2*pi/lambda;%波数
hole_ring_back=1e-3;
% hole_ring_back=0;
hole_ring_after=3e-3;
n=1;
R = n * 5 * 2 * lambda;%ROC曲率半径
a = n * 5 * lambda;%注意这里的a是孔径的一半
fnumber=R/(2*a);%所以f-number=曲率半径/孔径（2*a）
d = sqrt(R^2 - a^2);%理论焦点到孔径中心的距离
u=normal_velocity(P,R,a,0,medium.density,medium.soundspeed);

%划分网格点
xmin=-a;%观察点坐标的范围
xmax=-xmin;
ymax=0;
ymin=0;
zDiff=0.7*d;
zmin=R-zDiff;
zmax=R+zDiff;

dx = lambda/6; %网格点的步长
dz = lambda/6;

x=xmin:dx:xmax;%网格点的分布
z=zmin:dz:zmax;

nx=length(x);%网格点的点数
nz=length(z);

error_zeros_xz=zeros(nx,nz); %定义全零误差矩阵

%换能器离散化为点声源  注意离散化之后点声源的大小
dr=lambda/6;%球面分割成很多个圆环，对应的半径是r，dr是半径增加的步长
ndr=round((a-0)/dr);
dr=(a-0)/ndr;
r_back=0:dr:a-dr;%点声源前一段弧长对应的r
r_after=dr:dr:a;%点声源后一段弧长对应的r
r=r_after-dr/2;%第i个环带对应的中心点的r
Sm=lambda/6;%最外层环带（2pi*a）离散化后对应的dS的弧长小于lambda/6
median=fix(length(r)/2)+1;%fix截位取整
ntheta=round(2*pi*r(median)/Sm);%最外层一个环带的划分点数,取整
dtheta=2*pi./ntheta;%根据取整重新调整每个点声源的对应弧度
theta_after=dtheta:dtheta:2*pi;%每个环带离散成多个点对应的弧度数组
theta_back=0:dtheta:(2*pi-dtheta);%每个环带离散成多个点对应的弧度数组
theta=theta_after-dtheta/2;%第i个环带离散成点声源，dS中间的点对应的弧度数组
X=cos(theta)'*r;%点源的x坐标  矩阵形式 ntheta*nr
Y=sin(theta)'*r;%点源的y坐标
Z0=R-sqrt(R*R-r.*r);%点声源三维空间中的z坐标
Z=repmat(Z0,ntheta,1); % repmat( A , m , n )：将向量／矩阵在垂直方向复制m次，在水平方向复制n次。

%定义换能器法向阵速
u_ring_res=repmat(u,ntheta,ndr);
back_index=find(r_back<=hole_ring_back,1,'last');%判断环的内径落在哪个范围内
after_index=find(r_after>=hole_ring_after,1,'first');%判断环的外径落在哪个范围内
if back_index==after_index
    index_ring=back_index;
else index_ring=back_index:after_index;
end
u_ring_res(:,index_ring)=[0]; %将开孔部分法向阵速设为0

% 求y=0时xz平面
y=0;

%求开孔为环形剩余部分的声压
tic
for ix=1:nx  %观察网格点x方向的坐标
    for iz=1:nz  %观察网格点z方向的坐标
        rn=sqrt((X-x(ix)).^2+(Y-y).^2+(Z-z(iz)).^2);%观察点到点源的距离
        dS_ring=r.*dtheta.*R.*(asin(r_after./R)-asin(r_back./R));%每一环带离散的dS的大小，第i环离散的点声源面积dS=第i环的离散弧长*dr对应的短弧
        dS=repmat(dS_ring,ntheta,1); % repmat( A , m , n )：将向量／矩阵在垂直方向复制m次，在水平方向复制n次。
        A=dS.*exp(-1i.*k.*rn)./rn.*u_ring_res;
        B=sum(sum(A));%对上述求得的值累加
        pr_ring_res(ix,iz)=1i*medium.density*medium.soundspeed*k/(2*pi)*B; %乘以相关参数得到声压p
    end
end
toc

%求开孔部分的声压
%定义开孔部分法向阵速
u_ring_hole=zeros(ntheta,ndr);
u_ring_hole(:,index_ring)=[u]; %将开孔部分法向阵速设为0
tic
for ix=1:nx  %观察网格点x方向的坐标
    for iz=1:nz  %观察网格点z方向的坐标
        rn=sqrt((X-x(ix)).^2+(Y-y).^2+(Z-z(iz)).^2);%观察点到点源的距离
        dS_ring=r.*dtheta.*R.*(asin(r_after./R)-asin(r_back./R));%每一环带离散的dS的大小，第i环离散的点声源面积dS=第i环的离散弧长*dr对应的短弧
        dS=repmat(dS_ring,ntheta,1); % repmat( A , m , n )：将向量／矩阵在垂直方向复制m次，在水平方向复制n次。
        A=dS.*exp(-1i.*k.*rn)./rn.*u_ring_hole;
        B=sum(sum(A));%对上述求得的值累加
        pr_ring_hole(ix,iz)=1i*medium.density*medium.soundspeed*k/(2*pi)*B; %乘以相关参数得到声压p
    end
end
toc

%求完整换能器的声压
tic
for ix=1:nx  %观察网格点x方向的坐标
    for iz=1:nz  %观察网格点z方向的坐标
        rn=sqrt((X-x(ix)).^2+(Y-y).^2+(Z-z(iz)).^2);%观察点到点源的距离
        dS_ring=r.*dtheta.*R.*(asin(r_after./R)-asin(r_back./R));%每一环带离散的dS的大小，第i环离散的点声源面积dS=第i环的离散弧长*dr对应的短弧
        dS=repmat(dS_ring,ntheta,1); % repmat( A , m , n )：将向量／矩阵在垂直方向复制m次，在水平方向复制n次。
        A=dS.*exp(-1i.*k.*rn)./rn;
        B=sum(sum(A));%对上述求得的值累加
        pr_all(ix,iz)=1i*medium.density*u*medium.soundspeed*k/(2*pi)*B; %乘以相关参数得到声压p
    end
end
toc

%比较声压的幅值和相位
pr_all_amplitude=abs(pr_all);
pr_all_phase=angle(pr_all);
pr_sum=pr_ring_hole+pr_ring_res;
pr_sum_amplitude=abs(pr_sum);
pr_sum_phase=angle(pr_sum);
error_phase=abs(pr_sum_phase-pr_all_phase)./pr_all_phase;
error_amplitude=abs(pr_sum_amplitude-pr_all_amplitude)./pr_all_amplitude;
I_all=abs(pr_all).^2./(2*medium.soundspeed*medium.density);
I_hole=abs(pr_ring_hole).^2./(2*medium.soundspeed*medium.density);
I_res=abs(pr_ring_res).^2./(2*medium.soundspeed*medium.density);
max_hole=find_maxpoint(I_hole);
max_res=find_maxpoint(I_res);
max_all=find_maxpoint(I_all);


% %画图
figure(1);
surf(z*1000,x*1000,abs(pr_all));
axis equal;
shading interp;
colorbar;
xlabel('z (mm) ');
ylabel('x (mm) ');
title('the pressure amplitude of all(R=15mm,a=7.5mm,hole_a=1mm-3mm)');
% hold on;
figure(2);
surf(z*1000,x*1000,abs(pr_ring_res));
axis equal;
shading interp;
colorbar;
xlabel('z (mm) ');
ylabel('x (mm) ');
title('the pressure amplitude of residue(R=15mm,a=7.5mm,hole_a=1mm-3mm)');
figure(3);
surf(z*1000,x*1000,abs(pr_ring_hole));
axis equal;
shading interp;
colorbar;
xlabel('z (mm) ');
ylabel('x (mm) ');
title('the pressure amplitude of ring(R=15mm,a=7.5mm,hole_a=1mm-3mm)');
 figure(4);
surf(z*1000, x*1000,error_amplitude); 
axis equal;
shading interp;
colorbar;
xlabel('z (mm) ');
ylabel('x (mm) ');
title('the error of amplitude(R=15mm,a=7.5mm,hole_a=1mm-3mm)');
figure(5);
surf(z*1000, x*1000,error_phase);
axis equal;
shading interp;
colorbar;
xlabel('z (mm) ');
ylabel('x (mm) ');
title('the error of phase(R=15mm,a=7.5mm,hole_a=1mm-3mm))');
