clear all;
tic
load('pressure_xz_muscle1MPa.mat');
p_max=max(abs(p(:)));
[index]=find_maxpoint(abs(p));
xindex=index(1);
yindex=index(2);
zindex=index(3);
%  [focalsize_x,focalsize_y,focalsize_z] = focalsize(P3D,p_max,xi,yi,zi);%计算焦点大小

f=1e6;
% muscle;
% blood;
Alpha_muscle=9;
rho_muscle=1140;
c_muscle=1570;

XStep=dx;
deltaz=dz;

P3D=p;%(:,yindex,:);
Q3D = Alpha_muscle*abs(P3D).^2/(rho_muscle*c_muscle);

delta_t = 0.5;
k = 0.498;
C = 3480;
a1 = rho_muscle*C/delta_t + 2*k/(XStep*XStep)+k/(deltaz*deltaz);
b1 = -k/(2*XStep*XStep);
c1 = -k/(2*deltaz*deltaz);
a2 = rho_muscle*C/delta_t - 2*k/(XStep*XStep)-k/(deltaz*deltaz);
b2 = -b1;
c2 = -c1;
e = ones(size(P3D,1)*size(P3D,2)*size(P3D,3),1);

 Qc =0;% Wb * rho_b * Cb;

A = spdiags([c1*e,b1*e,b1*e,a1*e,b1*e,b1*e,c1*e],[-(nx*ny),-nx,-1,0,1,nx,(nx*ny)],length(e),length(e));

B = spdiags([c2*e,b2*e,b2*e,a2*e,b2*e,b2*e,c2*e],[-(nx*ny),-nx,-1,0,1,nx,(nx*ny)],length(e),length(e));

% Qb = Wb * rho_b * Cb * Tb;
% Qb = Qb*ones(size(P3D,1)*size(P3D,2)*size(P3D,3),1);

Q = reshape(Q3D,size(P3D,1)*size(P3D,2)*size(P3D,3),1);
% Q = Q + Qb;

T = 37*ones(size(P3D,1)*size(P3D,2)*size(P3D,3),1);
T_focus=[];
Tmax_near = [];
D = zeros(length(Q),1);
Q_perfusion = 0;    %Qc * T;
 
n = 0;
volume=0;

for n=1:10

[T,T_focus(n),D,Q_perfusion] = ThermalDose(A,B,D,T,Q,Q_perfusion,Qc,delta_t);%%%%%%热计量计算
     T_nearfield = reshape(T,size(P3D,1),size(P3D,2),size(P3D,3));%%%%%近场温升
    Tmax_near(n) = max(max(max(T_nearfield(:,:,1:60))));  %%%%%%%%%近场温升最大值
    
    TD = reshape(D,size(P3D,1),size(P3D,2),size(P3D,3));
    n = n + 1;
end
for i=1:length(D)
    if D(i)>240
        volume=volume+1;
    end
end

m = n;
time = delta_t*m;  

Trise = reshape(T,size(P3D,1),size(P3D,2),size(P3D,3));
Trise = Trise - 37*ones(size(P3D,1),size(P3D,2),size(P3D,3));

% Trise_xz=squeeze(Trise(:,yindex,:));
% Trise_xy=squeeze(Trise(:,:,zindex));

figure(1);
% surf(squeeze(Trise(:,:,zindex)));
surf(squeeze(Trise));
shading interp;
colormap(jet);
toc
% error=abs(Tnew-Trise_xy)./Tnew;  %和HIFU simulator做对比
% surf(error);
% shading interp;
% axis equal;
% colormap(jet);
