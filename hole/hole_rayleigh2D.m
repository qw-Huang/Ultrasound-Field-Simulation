%使用dS1点声源划分方法，在2d平面xz面，加开孔的换能器声强计算
clear all;
f0=1e6;%定义频率和声功率
P=100;
n=1;
medium = set_medium('lossless');%定义介质（单层：水），可以改成多层set_layered_medium
lambda = medium.soundspeed/f0;%波长=c/f
k=2*pi/lambda;%波数
hole_a=n*7.5e-4;
R = n * 5 * 2 * lambda;%ROC曲率半径
a = n * 5 * lambda;%注意这里的a是孔径的一半
fnumber=R/(2*a);%所以f-number=曲率半径/孔径（2*a）
d = sqrt(R^2 - a^2);%理论焦点到孔径中心的距离
u=normal_velocity(P,R,a,hole_a,medium.density,medium.soundspeed);

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
ndr=round((a-dr-hole_a)/dr);
dr=(a-dr-hole_a)/ndr;
r_back=(0+hole_a):dr:a-dr;%点声源前一段弧长对应的r
r_after=(dr+hole_a):dr:a;%点声源后一段弧长对应的r
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

% 求y=0时xz平面
y=0;

%rayleigh积分计算xz平面的声场  
tic
for ix=1:nx  %观察网格点x方向的坐标
    for iz=1:nz  %观察网格点z方向的坐标
        rn=sqrt((X-x(ix)).^2+(Y-y).^2+(Z-z(iz)).^2);%观察点到点源的距离
        dS_ring=r.*dtheta.*R.*(asin(r_after./R)-asin(r_back./R));%每一环带离散的dS的大小，第i环离散的点声源面积dS=第i环的离散弧长*dr对应的短弧
        dS=repmat(dS_ring,ntheta,1); % repmat( A , m , n )：将向量／矩阵在垂直方向复制m次，在水平方向复制n次。
        A=dS.*exp(-1i.*k.*rn)./rn;
        B=sum(sum(A));%对上述求得的值累加
        pr(ix,iz)=1i*medium.density*u*medium.soundspeed*k/(2*pi)*B; %乘以相关参数得到声压p
    end
end
toc

%比较dS的累加和解析面积
S_discrete=sum(sum(dS));%离散之后的全部面积累加
S_theory=2*pi*R*(R-d);%理论计算球面面积，即解析面积
dSi=sum(dS);%求每一环面积的累加，也就是把每一列（70个数）累加
di_back=sqrt(R^2-r_back.^2);%从0开始作为第i层环的d
di_after=sqrt(R^2-r_after.^2);%从dr开始作为第i层环的d
dSi_theory=2*pi*R*(R-di_after)-2*pi*R*(R-di_back);%第i环的环带面积
error_dSi=abs(dSi-dSi_theory)./dSi;%求解析解环带面积和数值解环带面积的误差

%声压转化为声场计算
I_pr=acousticintensity(pr,medium.density,medium.soundspeed); 

%画图
figure(1);
surf(z*1000,x*1000,I_pr); %横坐标为z，纵坐标为x，颜色代表p1的大小  
shading interp
title('Rayleigh(dS1 hole_a=7.5e-4)  ');
axis equal;%定义坐标的比例相同
colorbar
xlabel('z（mm） ');
ylabel('x (mm) ');
zlabel('normalized pressure');

% figure(2);
% surf(z*1000, x*1000, I_prs_nor);%对prs归一化   
% shading interp;
% colorbar;
% axis equal;
% title('Rayleigh(FOCUS)');
% xlabel('z (mm) ');
% ylabel('x (mm) ');
% zlabel('normalized pressure');

%对比两种方法画出来的结果
% figure(3);
% surf(z*1000, x*1000, error_zeros_xz); 
% shading interp %去掉网格，平滑曲面
% axis equal;
% colorbar  %加颜色条
% xlabel('z (mm) ');
% ylabel('x (mm) ');
% title('error between Rayleigh and Rayleigh (FOCUS)(R=150mm,a=75mm)');

figure(4);
histogram(dS_ring,10);
title('点源面积dS分布');

figure(5);
plot(error_dSi);
xlabel('ith ring ');
ylabel('error ');
title('error between the sum of dS and theory calculation');
