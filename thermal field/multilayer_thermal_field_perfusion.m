
 clear all;
 load('pressure3Dmultilayer_nohole.mat');
syms r;
f=f0;
delta_t=0.1;
delta_r=dx;   %1e-4;
delta_x=dz;  %1e-4;
% xx=x(98:417);   %只取组织部分 x =±0.04
xx=x(83:579);   %只取组织部分 x =±0.062
gridNum_x=length(z);  %1700;
gridNum_r=length(xx);  %1000;
times_heating=4;
nt_heating=times_heating/delta_t;
times_cooling=0;
nt_cooling=times_cooling/delta_t;
times=nt_heating+nt_cooling;
pressure_interest_area=squeeze(p_abs(83:579,y_index,:)); %只取组织部分 x =±0.04
pressure_xr=pressure_interest_area';%矩阵转置

Cb = 3800;  
T_ambient=310;
     
dBperNeper = 20 * log10(exp(1));
alpha_z2=medium2.attenuationdBcmMHz/dBperNeper*100*(f0/1e6)^1.1;
alpha_z3=medium3.attenuationdBcmMHz/dBperNeper*100*(f0/1e6)^1.1;
alpha_z4=medium4.attenuationdBcmMHz/dBperNeper*100*(f0/1e6)^1.1;
alpha_z5=medium5.attenuationdBcmMHz/dBperNeper*100*(f0/1e6)^1.1;
alpha_z2=repmat(alpha_z2,length(z2),gridNum_r);
alpha_z3=repmat(alpha_z3,length(z3),gridNum_r);
alpha_z4=repmat(alpha_z4,length(z4),gridNum_r);
alpha_z5=repmat(alpha_z5,length(z5),gridNum_r);
alpha=[alpha_z2;alpha_z3;alpha_z4;alpha_z5];

c_z2=repmat(medium2.soundspeed,length(z2),gridNum_r);
c_z3=repmat(medium3.soundspeed,length(z3),gridNum_r);
c_z4=repmat(medium4.soundspeed,length(z4),gridNum_r);
c_z5=repmat(medium5.soundspeed,length(z5),gridNum_r);
c=[c_z2;c_z3;c_z4;c_z5];

Kt_z2=repmat(medium2.thermalconductivity,length(z2),gridNum_r);
Kt_z3=repmat(medium3.thermalconductivity,length(z3),gridNum_r);
Kt_z4=repmat(medium4.thermalconductivity,length(z4),gridNum_r);
Kt_z5=repmat(medium5.thermalconductivity,length(z5),gridNum_r);
Kt=[Kt_z2;Kt_z3;Kt_z4;Kt_z5];

rou_z2=repmat(medium2.density,length(z2),gridNum_r);
rou_z3=repmat(medium3.density,length(z3),gridNum_r);
rou_z4=repmat(medium4.density,length(z4),gridNum_r);
rou_z5=repmat(medium5.density,length(z5),gridNum_r);
rou=[rou_z2;rou_z3;rou_z4;rou_z5];

Ct_z2=repmat(medium2.specificheat,length(z2),gridNum_r);
Ct_z3=repmat(medium3.specificheat,length(z3),gridNum_r);
Ct_z4=repmat(medium4.specificheat,length(z4),gridNum_r);
Ct_z5=repmat(medium5.specificheat,length(z5),gridNum_r);
Ct=[Ct_z2;Ct_z3;Ct_z4;Ct_z5];

Wb_z2=repmat(medium2.bloodperfusion,length(z2),gridNum_r);
Wb_z3=repmat(medium3.bloodperfusion,length(z3),gridNum_r);
Wb_z4=repmat(medium4.bloodperfusion,length(z4),gridNum_r);
Wb_z5=repmat(medium5.bloodperfusion,length(z5),gridNum_r);
Wb=[Wb_z2;Wb_z3;Wb_z4;Wb_z5];

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
% coefficients in BHTE
m = (2 * Kt * delta_t) ./ (rou .* Ct);
n = m / (12 * delta_x^2);
o = 2 * Wb .* Cb * delta_t ./ (rou .* Ct);
q = -5 / (2 * delta_r^2) - 5 / (2 * delta_x^2) - o ./ m; %血流灌注
s =2 * delta_t * alpha ./ ( Ct .* c .* (rou).^2);
% generate matrix Ta
Ta = T_ambient * ones(gridNum_x,gridNum_r);
T0 = T_ambient * ones(gridNum_x,gridNum_r);
T1 = T_ambient * ones(gridNum_x,gridNum_r);

input=ones(times,1);
if times_cooling>0
    input(nt_heating+1:times)=0;
end
D=zeros(gridNum_x*gridNum_r,1);
% calculate T 
for i=1:times
    Tp_2 = [T1(3:gridNum_x,:);zeros(2,gridNum_r)];
    Tp_1 = [T1(2:gridNum_x,:);zeros(1,gridNum_r)];
    Tm_1 = [zeros(1,gridNum_r);T1(1:gridNum_x - 1,:)];
    Tm_2 = [zeros(2,gridNum_r);T1(1:gridNum_x - 2,:)];

T2 =( m .* (T1 * A + q .* T1) + n .* (-Tp_2 + 16 * Tp_1 + 16 * Tm_1 - Tm_2)+o .* Ta  +4*T1 -T0 + input(i)*s .* abs(pressure_xr).^2)./3; %去掉血流灌注不考虑+o * Ta
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

       