% focus自带的rayleigh点数迭代对声场计算的影响
clc;clear all;
f0=1e6;u=1;
medium = set_medium('lossless');%改成set_layered_medium
lambda = medium.soundspeed/f0;
k=2*pi/lambda;
R = 5 * 2 * lambda;
a = 5 * lambda;
fnumber=R/(2*a);
d = sqrt(R^2 - a^2);
%划分网格点
xmin=-a;
xmax=a;
y=0;
ymax=0;
ymin=0;
zmin=R-d;
zmax=R+d;
nz = 101; % ok to sample the origin
nx = 61;
dx = 2 * a / (nx - 1);
dz = 2 * d / (nz - 1);
x=xmin:dx:xmax;
z=zmin:dz:zmax;

% Determine where the source pressure will be calculated
z0 = lambda/4;

 %调用focus内置
omega = 2 * pi * f0;
phi0 = asin(a/R);
xdcr = get_spherical_shell(a,R);

if nz > 1,
    dz = 2 * d / (nz - 1);
else
    dz = 0;
end
if nx > 1,
    dx = 2 * a / (nx - 1);
else
    dx = 0;
end

delta = [dx 0 dz];
ps = set_coordinate_grid(delta, xmin, xmax, ymin, ymax, zmin, zmax);
ndiv = 10:10:100;
dflag = 0;

tic
for i=1:length(ndiv)
    prs=rayleigh_cw(xdcr,ps,medium,ndiv(i),f0,dflag);
    nor_prs= abs(squeeze(prs))/max(max(abs(squeeze(prs))));
    prs2=rayleigh_cw(xdcr,ps,medium,ndiv(i)+100,f0,dflag);
    nor_prs2= abs(squeeze(prs2))/max(max(abs(squeeze(prs2))));
    error(i)=max(max(abs(prs2-prs)/abs(max(max(prs2)))));
end
toc

axis equal;
figure(1);
plot(ndiv,error);   
p2=abs(squeeze(prs));
xlabel('ndiv');
ylabel('normalized error');    
datestr(now)
