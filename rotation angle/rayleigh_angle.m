%换能器在xz平面偏转一定角度
clear all;
f0=1e6;%定义频率和法向阵速
P=100;
angle_rotation=pi/36;
medium = set_medium('muscle');%定义介质（单层：水），可以改成多层set_layered_medium
lambda = medium.soundspeed/f0;%波长=c/f
dBperNeper = 20 * log10(exp(1));
attenuationNeperspermeter=medium.attenuationdBcmMHz/dBperNeper*100*f0/1e6;
k=2*pi/lambda;%波数

R = 75e-3;%ROC曲率半径
a=30e-3;
fnumber=R/(2*a);%所以f-number=曲率半径/孔径（2*a）
d = sqrt(R^2 - a^2);%理论焦点到孔径中心的距离
u=normal_velocity(P,R,a,0,medium.density,medium.soundspeed);

%划分网格点
xmin=-1.5*a-10e-3;%观察点坐标的范围
xmax=1.5*a-10e-3;
ymax=0;
ymin=0;

zmin=10e-3;
% zmin=R-zDiff;
zmax=90e-3;

dx = lambda/6; %网格点的步长
dz = lambda/6;

x=xmin:dx:xmax;%网格点的分布
z=zmin:dz:zmax;

nx=length(x);%网格点的点数
nz=length(z);

pr_zeros_xz=zeros(nx,nz);
Xp=zeros(nx,nz);
Zp=zeros(nx,nz);

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
% X=cos(theta)'*r;%点源的x坐标  矩阵形式 ntheta*nr
X=cos(theta)'*r;
Y=sin(theta)'*r;%点源的y坐标
% Z0=R-sqrt(R*R-r.*r);%点声源三维空间中的z坐标
Z0=(R-sqrt(R*R-r.*r));
Z=repmat(Z0,ntheta,1); % repmat( A , m , n )：将向量／矩阵在垂直方向复制m次，在水平方向复制n次。
Rx=0;Rz=R; %绕几何焦点旋转，在xz平面内
Z_angle=cos(angle_rotation)*(Z-Rz)-sin(angle_rotation)*(X-Rx)+Rz; %二维空间中坐标旋转公式
X_angle=sin(angle_rotation)*(Z-Rz)+cos(angle_rotation)*(X-Rx)+Rx;

% 求y=0时xz平面
y=0;

%rayleigh积分计算xz平面的声场  
tic
for ix=1:nx  %观察网格点x方向的坐标
    for iz=1:nz  %观察网格点z方向的坐标
        rn=sqrt((X_angle-x(ix)).^2+(Y-y).^2+(Z_angle-z(iz)).^2);%观察点到点源的距离
        dS_ring=r.*dtheta.*R.*(asin(r_after./R)-asin(r_back./R));%每一环带离散的dS的大小，第i环离散的点声源面积dS=第i环的离散弧长*dr对应的短弧
        dS=repmat(dS_ring,ntheta,1); % repmat( A , m , n )：将向量／矩阵在垂直方向复制m次，在水平方向复制n次。
        A=dS.*exp(-1i.*k.*rn)./rn.*exp(-attenuationNeperspermeter.*rn);
        B=sum(sum(A));%对上述求得的值累加 
        pr(ix,iz)=1i*medium.density*u*medium.soundspeed*k/(2*pi)*B; %乘以相关参数得到声压p
    end
end
toc
pr_max=max(abs(pr(:)));
pr_abs=abs(pr);
%选择感兴趣区域
pr_nor=abs(pr)./pr_max;
[index]=find(pr_nor>=0.5);
pr_zeros_xz(index)=pr_nor(index);
X_all=repmat(x',1,nz);
Xp(index)=X_all(index);
Z_all=repmat(z,nx,1);
Zp(index)=Z_all(index);
%将感兴趣区域旋转回和z轴平行的位置，再进行比较
Zp_angle=cos(-angle_rotation)*(Zp-Rz)-sin(-angle_rotation)*(Xp-Rx)+Rz; 
Xp_angle=sin(-angle_rotation)*(Zp-Rz)+cos(-angle_rotation)*(Xp-Rx)+Rx;
Xp(index)=Xp_angle(index);%坐标为将焦区转回长轴和z轴重合
Zp(index)=Zp_angle(index);
axial=max(Zp(:))-min(Zp(:));%有问题 最小是0
radial=max(Xp(:))-min(Xp(:));

figure(3)
surf(z*1000,x*1000,abs(pr)./pr_max);
shading interp;
colorbar;
axis equal;
