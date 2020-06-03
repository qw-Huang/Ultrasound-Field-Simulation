function [temperature] = Thermal_field(dx,dz,dt,medium,times,pressure_xz)
%FDTD_THERMAL_FIELD 热场温度计算函数 单介质，不考虑血流灌注,输出温度
%  Thermal_field(dx,dz,dt,medium,times,pressure_xz)

pressure1=pressure_xz'; %input pressure
syms r;
delta_r=dx;   %1e-4;
delta_t=dt;
delta_x=dz;  %1e-4;
gridNum_x=length(z);  %1700;
gridNum_r=nx;  %1000;
%  Convert attenuation to Np/m
dBperNeper = 20 * log10(exp(1));
attenuationNeperspermeter=medium.attenuationdBcmMHz/dBperNeper*100*(f0/1e6)^1.1;
alpha=attenuationNeperspermeter;
Cb = 3800;  
T_ambient=310;
%% BHTE parameters
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
% coefficients in BHTE
m = (2 * medium.Kt * delta_t) / (medium.rou * medium.Ct);
n = m / (12 * delta_x^2);
% o = 2 * medium.Wb * Cb * delta_t / (medium.rou * medium.Ct);
q = -5 / (2 * delta_r^2) - 5 / (2 * delta_x^2);% - o / m去掉血流灌注
s =2 * delta_t * alpha / ( medium.Ct * medium.c * (medium.rou)^2);
% generate matrix Ta
% Ta = T_ambient * ones(gridNum_x,gridNum_r);
T0 = T_ambient * ones(gridNum_x,gridNum_r);
T1 = T_ambient * ones(gridNum_x,gridNum_r);

%% calculate T 
for i=1:times
    Tp_2 = [T1(3:gridNum_x,:);zeros(2,gridNum_r)];
    Tp_1 = [T1(2:gridNum_x,:);zeros(1,gridNum_r)];
    Tm_1 = [zeros(1,gridNum_r);T1(1:gridNum_x - 1,:)];
    Tm_2 = [zeros(2,gridNum_r);T1(1:gridNum_x - 2,:)];
T2 =( m * (T1 * A + q * T1) + n * (-Tp_2 + 16 * Tp_1 + 16 * Tm_1 - Tm_2)  +4*T1 -T0 + s * abs(pressure1).^2)./3; %去掉血流灌注不考虑+o * Ta
% def clearBoundary(T2):
    T2(:,1) = T_ambient * ones(gridNum_x,1);
    T2(:,2) = T_ambient * ones(gridNum_x,1);
    T2(:,gridNum_r - 1) = T_ambient * ones(gridNum_x,1);
    T2(:,gridNum_r)     = T_ambient * ones(gridNum_x,1);
    T2(1,:) = T_ambient * ones(1,gridNum_r);
    T2(2,:) = T_ambient * ones(1,gridNum_r);
    T2(gridNum_x - 1,:)= T_ambient * ones(1,gridNum_r);
    T2(gridNum_x ,:) = T_ambient * ones(1,gridNum_r);    
%     T_focus(i)=T2(140,127); calculate focus temeprature changes
    T0=T1;
    T1=T2;  
end
%% 
temperature = T2;     %output T2 gridNum_x*gridNum_r
end

