%使用内置rayleigh计算xy面声场分布，结果存在问题xy平面
clc;
clear all;
clear all;
f0=1e6;%定义频率和法向阵速
u=1;
medium = set_medium('water');%定义介质（单层：水），可以改成多层set_layered_medium
lambda = medium.soundspeed/f0;%波长=c/f

k=2*pi/lambda;%波数

R = 75e-3;%ROC曲率半径
a = 30e-3;%注意这里的a是孔径的一半
fnumber=R/(2*a);%所以f-number=曲率半径/孔径（2*a）
d = sqrt(R^2 - a^2);%理论焦点到孔径中心的距离

%划分网格点
xmin=-1.5*a;%观察点坐标的范围
xmax=-xmin;
ymin=xmin;
ymax=-ymin;
z0=30e-3;

dx = lambda/6; %网格点的步长
dy = lambda/6; 
dz = 0;

x=xmin:dx:xmax;%网格点的分布
y=ymin:dy:ymax;
z=z0;

nx=length(x);%网格点的点数
ny=length(y);
% error_zeros_xz=zeros(nx,nz); %定义全零误差矩阵

% %调用focus内置rayleigh函数rayleigh_cw()
% omega = 2 * pi * f0;
% phi0 = asin(a/R);
% xdcr=get_spherical_shell(a,R);
% 
% delta = [dx dy 0];%网格点步长
% ps = set_coordinate_grid(delta, xmin, xmax, ymin, ymax, z0, z0);%网格点的划分
% ndiv = 100;%积分的点数
% dflag = 0;
% 
% tic
% prs=fnm_call(xdcr,ps,medium,ndiv,f0);%FOCUS自带的rayleigh计算
% toc
% I_prs=acousticintensity(prs,medium.density,medium.soundspeed);
% I_prs_nor=I_prs./max(I_prs(:));

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
Z0=R-sqrt(R*R-r.*r);%点声源三维空间中的z坐标
Z=repmat(Z0,ntheta,1); % repmat( A , m , n )：将向量／矩阵在垂直方向复制m次，在水平方向复制n次。

%rayleigh积分计算xz平面的声场  
tic
for ix=1:nx  %观察网格点x方向的坐标
    for iy=1:ny  %观察网格点z方向的坐标
        rn=sqrt((X-x(ix)).^2+(Y-y(iy)).^2+(Z-z0).^2);%观察点到点源的距离
        dS_ring=r.*dtheta.*R.*(asin(r_after./R)-asin(r_back./R));%每一环带离散的dS的大小，第i环离散的点声源面积dS=第i环的离散弧长*dr对应的短弧
        dS=repmat(dS_ring,ntheta,1); % repmat( A , m , n )：将向量／矩阵在垂直方向复制m次，在水平方向复制n次。
        A=dS.*exp(-1i.*k.*rn)./rn;
        B=sum(sum(A));%对上述求得的值累加
        pr(ix,iy)=1i*medium.density*u*medium.soundspeed*k/(2*pi)*B; %乘以相关参数得到声压p
    end
end
toc
% 
% I_pr=acousticintensity(pr,medium.density,medium.soundspeed); 
% I_pr_nor=I_pr./max(I_pr(:));
% error=abs(I_pr_nor-I_prs_nor)./I_pr_nor;

% % plot the pressures
% figure(1);
% surf(y*1000,x*1000,I_prs_nor);%对prs归一化   
% shading interp
% colorbar;
% axis equal;
% xlabel('y(mm)');
% ylabel('x(mm)');
% 
% 
% figure(2);
% surf(y*1000,x*1000,I_pr_nor);%对pr归一化   
% shading interp;
% axis equal;
% colorbar;
% % title('z=R,Rayleigh Sommerfeld Result');
% xlabel('y(mm)');
% ylabel('x(mm)');
% 
% figure(3);
% surf(y*1000,x*1000,error);%对pr归一化   
% shading interp;
% axis equal;
% colorbar;
% xlabel('y(mm)');
% ylabel('x(mm)');