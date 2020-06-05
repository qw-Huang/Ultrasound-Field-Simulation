clear all;
load('pressure3Dmuscle.mat');
syms r;
% f=f0;
delta_t=0.1;   %(unit:s)
delta_r=dx;   %1e-4;
delta_x=dz;  %1e-4;
xx=x(98:417);   %只取组织部分 x =±0.04
gridNum_x=length(z);  %1700;
gridNum_r=length(xx);  %1000;

pressure_interest_area=squeeze(p_abs(98:417,y_index,:)); %只取组织部分 x =±0.04
pressure_xr=pressure_interest_area';%矩阵转置

times=20;
Cb = 3800;  
T_ambient=310;

tic
c1 = -1 / (12 * delta_r^2);
d1 = 4 / (3 * delta_r^2);
e1 = d1;
f1 = c1;

  c2(r)= -1 / (12 * (r + 0.5) * delta_r^2);
  
  d2(r)= 2 / (3 * (r + 0.5) * delta_r^2);

  e2(r)= -2 / (3 * (r + 0.5) * delta_r^2);
    
  f2(r)= 1 / (12 * (r + 0.5) * delta_r^2);
% generate matrix A
A = zeros(gridNum_r,gridNum_r);
for i =1:gridNum_r
        rj=abs(i-round(gridNum_r/2));
    if(i >= 2)
        A(i-1,i) = e1 + e2(rj);
    end
    if(i >= 3)
        A(i-2,i) = f1 + f2(rj);
    end
    if(i <= gridNum_r - 1)
        A(i+1,i) = d1 + d2(rj);
    end
    if(i <= gridNum_r - 2)
        A(i+2,i) = c1 + c2(rj);
    end
end

Kt=medium2.thermalconductivity;
rou=medium2.density;
Ct=medium2.specificheat;
c=medium2.soundspeed;
dBperNeper = 20 * log10(exp(1));
alpha=medium2.attenuationdBcmMHz/dBperNeper*100*(f0/1e6)^1.1;
% alpha=attenuationNeperspermeter;

% coefficients in BHTE
m = (2 * Kt * delta_t) / (rou * Ct);
n = m / (12 * delta_x^2);
% o = 2 * medium.Wb * Cb * delta_t / (medium.rou * medium.Ct);
q = -5 / (2 * delta_r^2) - 5 / (2 * delta_x^2); %- o / m; %血流灌注
s =2 * delta_t * alpha / (Ct * c * (rou)^2);
%  generate matrix Ta
% Ta = T_ambient * ones(gridNum_x,gridNum_r);
T0 = T_ambient * ones(gridNum_x,gridNum_r);
T1 = T_ambient * ones(gridNum_x,gridNum_r);

% pressure_xr=pressure_xz';
D=zeros(gridNum_x*gridNum_r,1);

% calculate T 
for i=1:times
    Tp_2 = [T1(3:gridNum_x,:);zeros(2,gridNum_r)];
    Tp_1 = [T1(2:gridNum_x,:);zeros(1,gridNum_r)];
    Tm_1 = [zeros(1,gridNum_r);T1(1:gridNum_x - 1,:)];
    Tm_2 = [zeros(2,gridNum_r);T1(1:gridNum_x - 2,:)];

T2 =( m * (T1 * A + q * T1) + n * (-Tp_2 + 16 * Tp_1 + 16 * Tm_1 - Tm_2)  +4*T1 -T0 + s * abs(pressure_xr).^2)./3; %去掉血流灌注不考虑+o * Ta
% def clearBoundary(T2):
    T2(:,1) = T_ambient * ones(gridNum_x,1);
    T2(:,2) = T_ambient * ones(gridNum_x,1);
    T2(:,gridNum_r - 1) = T_ambient * ones(gridNum_x,1);
    T2(:,gridNum_r)     = T_ambient * ones(gridNum_x,1);
    T2(1,:) = T_ambient * ones(1,gridNum_r);
    T2(2,:) = T_ambient * ones(1,gridNum_r);
    T2(gridNum_x - 1,:)= T_ambient * ones(1,gridNum_r);
    T2(gridNum_x ,:) = T_ambient * ones(1,gridNum_r);
    
    T_focus(i)=T2(z_index,round(gridNum_r/2));
    T0=T1;
    T1=T2;
%calculate thermal dose
T2re=reshape(T2,gridNum_x*gridNum_r,1); %二维数组重构为1维
T=T2re-273;
    for j=1:length(T)
        if T(j)>=43
            D(j)=D(j)+0.5^(43-T(j))*delta_t;
        else
            D(j)=D(j)+0.25^(43-T(j))*delta_t;
        end
    end

end
ThermalDose=reshape(D,gridNum_x,gridNum_r);
toc
figure(1);
surf(xx,z,pressure_xr);
shading interp;
axis equal;
figure ()
surf(T2-273);
shading interp;
axis equal;
colormap(jet);
figure();
plot(T_focus);
figure();
contourf(ThermalDose',[240,240]);

       