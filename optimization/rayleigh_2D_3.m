%测试计算时间 对代码进行改进，将两个 for循环改为1个for循环

clear all;
profile on;

f0=1e6;u=1;%定义频率和法向阵速
medium = set_medium('lossless');%定义介质（单层：水），可以改成多层set_layered_medium
lambda = medium.soundspeed/f0;%波长=c/f
k=2*pi/lambda;%波数

R = 4 * 5 * 2 * lambda;%ROC曲率半径
a = 4 * 5 * lambda;%注意这里的a是孔径的一半
fnumber=R/(2*a);%所以f-number=曲率半径/孔径（2*a）
d = sqrt(R^2 - a^2);%理论焦点到孔径中心的距离

%划分网格点
xmin=-a;%观察点坐标的范围
xmax=0;
ymax=-a;
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

%换能器离散化为点声源  注意离散化之后点声源的大小
dr=lambda/6;%球面分割成很多个圆环，对应的半径是r，dr是半径增加的步长
nr=round(a/dr);%球面离散的环带个数
dr=a/nr;%取整重新调整
r_back=0:dr:a-dr;%点声源前一段弧长对应的r
r_after=dr:dr:a;%点声源后一段弧长对应的r
r=r_after-dr/2;%第i个环带对应的中心点的r
Sm=lambda/6;%最外层环带（2pi*a）离散化后对应的dS的弧长小于lambda/6
median=length(r)/2+1; %取中间环带的索引
ntheta=round(2*pi*r(median)/Sm);%中间层一个环带的划分点数,取整
dtheta=2*pi./ntheta;%根据取整重新调整每个点声源的对应弧度
theta_after=dtheta:dtheta:2*pi;%每个环带离散成多个点对应的弧度数组
theta_back=0:dtheta:(2*pi-dtheta);%每个环带离散成多个点对应的弧度数组
theta=theta_after-dtheta/2;%第i个环带离散成点声源，dS中间的点对应的弧度数组
X_source=cos(theta)'*r;%点源的x坐标  矩阵形式 ntheta*nr
X_row=X_source(:)';%点声源在三维空间的横坐标，每行代表环数，每列代表一环离散的点数，按列展开展成一行，一行是按从内到外的环数排开
Y_source=sin(theta)'*r;%点源的y坐标
Y_row=Y_source(:)';%每行代表环数，每列代表一环离散的点数，按列展开展成一行，一行是按从内到外的环数排开
Z0=R-sqrt(R*R-r.*r);%点声源三维空间中的z坐标
Z_source=repmat(Z0,ntheta,1); % repmat( A , m , n )：将向量／矩阵在垂直方向复制m次，在水平方向复制n次。
Z_row=Z_source(:)';
% 求y=0时xz平面
y=0;

%把X Y Z分别扩展
X=repmat(X_row,nx,1);
Y=repmat(Y_row,nx,1);
Z=repmat(Z_row,nx,1);

x_colnum=x';%把x方向上网格点的坐标变成列向量
x_repeat=repmat(x_colnum,1,ntheta*nr);%将列向量扩展和X Y Z相同的矩阵
%rayleigh积分计算xz平面的声场  
tic

    for iz=1:nz  %观察网格点z方向的坐标
        rn=sqrt((X-x_repeat).^2+(Y-y).^2+(Z-z(iz)).^2);%观察点到点源的距离
        dS_ring=r.*dtheta.*R.*(asin(r_after./R)-asin(r_back./R));%每一环带离散的dS的大小，第i环离散的点声源面积dS=第i环的离散弧长*dr对应的短弧
        dS_source=repmat(dS_ring,ntheta,1); % repmat( A , m , n )：将向量／矩阵在垂直方向复制m次，在水平方向复制n次。
        dS_row=dS_source(:)';%把dS按照上面的方式展成一行
        dS=repmat(dS_row,nx,1);%dS扩展为和XYZ相同的矩阵
        A=dS.*exp(-1i*k.*rn)./rn;
        B=sum(A,2);%对上述求得的值累加
        pr_row=1i*medium.density*u*medium.soundspeed*k/(2*pi)*B; %乘以相关参数得到声压p
        pr_up(:,iz)=pr_row;
    end

toc
pr_down=flipud(pr_up);%对声压矩阵进行反转，相当于关于x=0轴对称
pr_all=[pr_up;pr_down];%将上述两个矩阵拼接起来，但是中间一行会复制两遍
pr=pr_all;
[row,column]=size(pr);
index_median=row/2;%取中间行
pr(index_median,:)=[]; %删掉中间重复的一行

%声压转化为声场计算
I_pr=acousticintensity(pr,medium.density,medium.soundspeed); 
I_pr_nor=I_pr./max(I_pr(:));
max_index=find_maxpoint(I_pr);%返回最大点位置坐标
profile viewer