%测试时间 对三维的代码进行改进，只计算四分之一个象限 
close all;
clear all;
clear all;
profile on;
f0=1e6;%定义频率和法向阵速
P=100;
medium = set_medium('lossless');%定义介质：单层->水，可改成多层，用函数set_layered_medium
lambda = medium.soundspeed/f0;%波长=c/f
k=2*pi/lambda;%波数

R =10* 5 * 2 * lambda;%ROC曲率半径
a =10* 5 * lambda;%注意这里a是孔径的一半
fnumber=R/(2*a);%所以f-number=曲率半径/孔径（2*a）
d = sqrt(R^2 - a^2);%理论焦点到孔径中心的距离
u=normal_velocity(P,R,a,0,medium.density,medium.soundspeed);

%划分网格点
xmin=-0.005;%观察点坐标的范围
xmax=0;
ymin=-0.005;
ymax=0;
zDiff=0.01;
zmin=R-zDiff;
zmax=R+zDiff;

dx = lambda/6; %网格点的步长
dy = lambda/6; 
dz = lambda/6;

x=xmin:dx:xmax;%网格点的分布
y=ymin:dy:ymax;
z=zmin:dz:zmax;

nx=length(x);%网格点的点数
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

%rayleigh积分计算三维空间的声场  
tic
for ix=1:nx  %观察网格点x方向的坐标划分
    for iy=1:ny %观察网格点y方向的坐标划分
       for iz=1:nz  %观察网格点z方向的坐标划分
            rn=sqrt((X-x(ix)).^2+(Y-y(iy)).^2+(Z-z(iz)).^2);%观察点到点源的距离
            dS_ring=r.*dtheta.*R.*(asin(r_after./R)-asin(r_back./R));%每一环带离散的dS的大小，第i环离散的点声源面积dS=第i环的离散弧长*dr对应的短弧
            dS=repmat(dS_ring,ntheta,1); % repmat( A , m , n )：将向量／矩阵在垂直方向复制m次，在水平方向复制n次。
            A=dS.*exp(-1i.*k.*rn)./rn;
            B=sum(sum(A));%对上述求得的值累加
            pr_1(ix,iy,iz)=1i*medium.density*u*medium.soundspeed*k/(2*pi)*B; %乘以相关参数得到声压p   
        end
    end
end
toc
pr_3=flipud(pr_1);%对声压矩阵进行反转，相当于关于x=0轴对称
pr_13=[pr_1;pr_3];%将上述两个矩阵拼接起来，但是中间一行会复制两遍
[row,column]=size(pr_13);
row_median=row/2;%取中间行
pr_13(row_median,:,:)=[]; %删掉中间重复的一行
pr_24=fliplr(pr_13);%对声压矩阵进行左右反转，相当于关于y=0轴对称
pr_all=[pr_13,pr_24];%将上述两个矩阵拼接起来，但是中间一行会复制两遍
[row,column,c]=size(pr_all);
col_median=column/2;%取中间列
pr_all(:,col_median,:)=[]; %删掉中间重复的一行

%声压转化为声场计算
I_pr=acousticintensity(pr_all,medium.density,medium.soundspeed); 
I_pr_nor=I_pr./max(I_pr(:));
max_index=find_maxpoint(I_pr);%返回最大点位置坐标
y_index=max_index(2);
x=-0.005:dx:0.005;
surf(z*1000,x*1000,squeeze(I_pr(:,y_index,:)));
shading interp;
axis equal;
profile viewer



