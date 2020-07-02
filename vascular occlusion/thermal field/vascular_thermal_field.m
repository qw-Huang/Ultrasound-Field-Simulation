% clear all;
% load('pressure3Dmulti_blood.mat');
syms r;
f=f0;
delta_t=1e-2;
delta_tb=1e-2;
Nt=round(delta_t/delta_tb);
delta_r=dx;   %1e-4;
delta_x=dz;  %1e-4;
% xx=x(98:417);   %只取组织部分 x =±0.04
xx=x(83:579);   %只取组织部分 x =±0.062
gridNum_x=length(z);  %1700;
gridNum_r=length(xx);  %1000;
times_heating=5;
nt_heating=round(times_heating/delta_t);
times_cooling=0;
nt_cooling=times_cooling/delta_t;
times=nt_heating+nt_cooling;
pressure_interest_area=squeeze(p_abs(83:579,y_index,:)); %只取组织部分 x =±0.04
pressure_xr=pressure_interest_area';%矩阵转置

% t_c = 0;	% cool-off duration (s)
% 
% % % periodic sonications
% n_c = 1;		% number of pulse cycles 重复脉冲数
% D = 50;		% duty cycle (%) 占空比
% t_p = 1;	% pulse cycle period (s) 脉冲循环周期
% 
% % determine number of integration steps at each stage:
% % N_i = round(t_i/delta_t);
% if(n_c==0)
%   N_p=0;
% else
%   N_p = round(t_p/delta_t);
% end
% N_c = round(t_c/delta_t);
% N = n_c*N_p+N_c;		% total number of integration steps
% Ts = delta_t*N;			% total simulation duration
% t=0:delta_t:Ts;
% % build input vector, which catalogs the HIFU beam on/off operation 
% if(N_p~=0)
%     pulse = zeros(1,N_p);
%   for n=1:round((1-0.01*D)*N_p)
%     pulse(n) = 1;
%   end
%   pulses = pulse;
%   for m=2:n_c
%     pulses = [pulses,pulse];
%   end
% end
% if(N_c~=0)
%   cooloff = zeros(1,N_c);
% end 
% % if(N_i~=0)
% %   initial = ones(1,N_i);
% % end
% 
% for n=1:n_c*N_p
%   input(n) = pulses(n);
% end
% for n=n_c*N_p+1:N
%     input(n) = cooloff(n-n_c*N_p);
% end
% 
% times=N;
Cb = 3800;  
T_ambient=310;
r0=0.0025;   % radius of vascular =2.5mm
V0=0.0067;   % average blood flow velocity 30cm/s
theta=pi/2;

dBperNeper = 20 * log10(exp(1));
alpha_z2=medium2.attenuationdBcmMHz/dBperNeper*100*(f0/1e6)^1.1;
alpha_z3=medium3.attenuationdBcmMHz/dBperNeper*100*(f0/1e6)^1.1;
alpha_z4=medium4.attenuationdBcmMHz/dBperNeper*100*(f0/1e6)^1.1;
alpha_z5=medium5.attenuationdBcmMHz/dBperNeper*100*(f0/1e6)^1.1;
alpha_z6=medium6.attenuationdBcmMHz/dBperNeper*100*(f0/1e6)^1.1;
alpha_z2=repmat(alpha_z2,length(z2),gridNum_r);
alpha_z3=repmat(alpha_z3,length(z3),gridNum_r);
alpha_z4=repmat(alpha_z4,length(z4),gridNum_r);
alpha_multi=[alpha_z2;alpha_z3;alpha_z4];
alpha_blood=alpha_z5;
alpha_muscle=alpha_z6;

c_z2=repmat(medium2.soundspeed,length(z2),gridNum_r);
c_z3=repmat(medium3.soundspeed,length(z3),gridNum_r);
c_z4=repmat(medium4.soundspeed,length(z4),gridNum_r);
c_multi=[c_z2;c_z3;c_z4];
c_blood=medium5.soundspeed;
c_muscle=medium6.soundspeed;

Kt_z2=repmat(medium2.thermalconductivity,length(z2),gridNum_r);
Kt_z3=repmat(medium3.thermalconductivity,length(z3),gridNum_r);
Kt_z4=repmat(medium4.thermalconductivity,length(z4),gridNum_r);
Kt_multi=[Kt_z2;Kt_z3;Kt_z4];
Kt_blood=medium5.thermalconductivity;
Kt_muscle=medium6.thermalconductivity;

rou_z2=repmat(medium2.density,length(z2),gridNum_r);
rou_z3=repmat(medium3.density,length(z3),gridNum_r);
rou_z4=repmat(medium4.density,length(z4),gridNum_r);
rou_multi=[rou_z2;rou_z3;rou_z4];
rou_blood=medium5.density;
rou_muscle=medium6.density;

Ct_z2=repmat(medium2.specificheat,length(z2),gridNum_r);
Ct_z3=repmat(medium3.specificheat,length(z3),gridNum_r);
Ct_z4=repmat(medium4.specificheat,length(z4),gridNum_r);
Ct_multi=[Ct_z2;Ct_z3;Ct_z4];
Ct_blood=medium5.specificheat;
Ct_muscle=medium6.specificheat;

Wb_z2=repmat(medium2.bloodperfusion,length(z2),gridNum_r);
Wb_z3=repmat(medium3.bloodperfusion,length(z3),gridNum_r);
Wb_z4=repmat(medium4.bloodperfusion,length(z4),gridNum_r);
Wb_multi=[Wb_z2;Wb_z3;Wb_z4];
Wb_blood=medium5.bloodperfusion;
Wb_muscle=medium6.bloodperfusion;

%index of interface
index_tissueblood=size(Ct_multi,1);
index_bloodmuscle=index_tissueblood+length(z5);

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
% coefficients in BHTE  multilayer
m_multi = (2 * Kt_multi * delta_t) ./ (rou_multi .* Ct_multi);
n_multi = m_multi / (12 * delta_x^2);
o_multi = 2 * Wb_multi .* Cb * delta_t ./ (rou_multi .* Ct_multi);
q_multi = -5 / (2 * delta_r^2) - 5 / (2 * delta_x^2) - o_multi ./ m_multi; %血流灌注
s_multi =2 * delta_t * alpha_multi ./ ( Ct_multi .* c_multi .* (rou_multi).^2);
%coefficients in BHTE blood
m_blood = 2 * Kt_blood * delta_tb / (rou_blood * Ct_blood);
n_blood = m_blood / (12 * delta_x^2);
q_blood = -5 / (2 * delta_r^2) - 5 / (2 * delta_x^2) ; %血流灌注
s_blood =2 * delta_tb * alpha_blood / ( Ct_blood * c_blood * (rou_blood) ^2);
%coefficients in BHTE muscle
m_muscle = (2 * Kt_muscle * delta_t) / (rou_muscle * Ct_muscle);
n_muscle = m_muscle / (12 * delta_x^2);
o_muscle = 2 * Wb_muscle * Cb * delta_t / (rou_muscle * Ct_muscle);
q_muscle = -5 / (2 * delta_r^2) - 5 / (2 * delta_x^2) - o_muscle / m_muscle; %血流灌注
s_muscle =2 * delta_t * alpha_muscle / ( Ct_muscle * c_muscle * (rou_muscle) ^2);

% generate matrix Ta
Ta = T_ambient * ones(gridNum_x,gridNum_r);
T0 = T_ambient * ones(gridNum_x,gridNum_r);
T1 = T_ambient * ones(gridNum_x,gridNum_r);

% generate different medium Ta and pressure
pressure_xr_multi=pressure_xr(1:index_tissueblood,:);
pressure_xr_blood=pressure_xr(index_tissueblood:index_bloodmuscle,:);
pressure_xr_muscle=pressure_xr(index_bloodmuscle:gridNum_x,:);
Ta_multi=Ta(1:index_tissueblood,:);
Ta_blood=Ta(index_tissueblood:index_bloodmuscle,:);
Ta_muscle=Ta(index_bloodmuscle:gridNum_x,:);

% radius of vascular  
index_vascular_center=round((nz5-1)/2);
z_blood=z(index_tissueblood:index_bloodmuscle);
rb=abs(z_blood-z5(index_vascular_center));
rb=repmat(rb',1,gridNum_r);

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
%calculate mulilayer thermal field
    T0_multi=T0(1:index_tissueblood,:);
    T1_multi=T1(1:index_tissueblood,:);
    Tp_2_multi=Tp_2(1:index_tissueblood,:);
    Tp_1_multi=Tp_1(1:index_tissueblood,:);
    Tm_2_multi=Tm_2(1:index_tissueblood,:);
    Tm_1_multi=Tm_1(1:index_tissueblood,:);
    
    T2_multi =( m_multi .* (T1_multi * A + q_multi .* T1_multi) + n_multi .* (-Tp_2_multi + 16 * Tp_1_multi + 16 * Tm_1_multi - Tm_2_multi)...
            +o_multi .* Ta_multi +4*T1_multi -T0_multi + input(i)*s_multi .* abs(pressure_xr_multi).^2)./3; %去掉血流灌注不考虑+o * Ta
%calculate blood thermal field
    T0_blood=T0(index_tissueblood:index_bloodmuscle,:);
    T1_blood=T1(index_tissueblood:index_bloodmuscle,:);
%     for j=1:Nt
    Tp_2_blood=Tp_2(index_tissueblood:index_bloodmuscle,:);
    Tp_1_blood=Tp_1(index_tissueblood:index_bloodmuscle,:);
    Tm_2_blood=Tm_2(index_tissueblood:index_bloodmuscle,:);
    Tm_1_blood=Tm_1(index_tissueblood:index_bloodmuscle,:);
%     nx_blood=size(T1_blood,1);
%     nx_multi=size(T1_multi,1);
%     Tp_2_blood = [T1_blood(3:nx_blood,:);zeros(2,gridNum_r)];
%     Tp_1_blood = [T1_blood(2:nx_blood,:);zeros(1,gridNum_r)];
%     Tm_1_blood = [zeros(1,gridNum_r);T1_blood(1:nx_blood - 1,:)];
%     Tm_2_blood = [zeros(2,gridNum_r);T1_blood(1:nx_blood - 2,:)];    

    Tq_2=[T1_blood(:,3:gridNum_r) zeros(size(T1_blood,1),2)];
    Tq_1=[T1_blood(:,2:gridNum_r) zeros(size(T1_blood,1),1)];
    Tn_1=[zeros(size(T1_blood,1),1) T1_blood(:,1:gridNum_r - 1)];
    Tn_2=[zeros(size(T1_blood,1),2) T1_blood(:,1:gridNum_r - 2)];

    T2_blood =( m_blood * (T1_blood * A + q_blood * T1_blood) + n_blood * (-Tp_2_blood + 16 * Tp_1_blood + 16 * Tm_1_blood - Tm_2_blood)...
            -4*delta_tb*V0*(1-(rb./r0).^2).*(sin(theta)*1/(12*delta_r)*(-Tq_2 + 8 * Tq_1 - 8 * Tn_1 + Tn_2)...
            +cos(theta)*1/(12*delta_x)*(-Tp_2_blood + 8 * Tp_1_blood - 8 * Tm_1_blood + Tm_2_blood)) ...
            +4*T1_blood -T0_blood + input(i)* s_blood * abs(pressure_xr_blood).^2)./3; 
%              +4*T1_blood -T0_blood + s_blood * abs(pressure_xr_blood).^2)./3; 
% % def clearBoundary(T2_blood):
%     T2_blood(:,1) = T_ambient * ones(nx_blood,1);
%     T2_blood(:,2) = T_ambient * ones(nx_blood,1);
%     T2_blood(:,gridNum_r - 1) = T_ambient * ones(nx_blood,1);
%     T2_blood(:,gridNum_r)     = T_ambient * ones(nx_blood,1);
%     T2_blood(1,:) = T2_multi(nx_multi,gridNum_r);
% %     T2_blood(gridNum_x ,:) = T2(index_bloodmuscle,gridNum_r);
% 
%         T0_blood=T1_blood;
%         T1_blood=T2_blood;
%     end
%calculate muscle thermal field
    T0_muscle=T0(index_bloodmuscle:gridNum_x,:);
    T1_muscle=T1(index_bloodmuscle:gridNum_x,:);
    Tp_2_muscle=Tp_2(index_bloodmuscle:gridNum_x,:);
    Tp_1_muscle=Tp_1(index_bloodmuscle:gridNum_x,:);
    Tm_2_muscle=Tm_2(index_bloodmuscle:gridNum_x,:);
    Tm_1_muscle=Tm_1(index_bloodmuscle:gridNum_x,:);
    
    T2_muscle =( m_muscle * (T1_muscle * A + q_muscle * T1_muscle) + n_muscle * (-Tp_2_muscle + 16 * Tp_1_muscle + 16 * Tm_1_muscle - Tm_2_muscle)...
            +o_muscle * Ta_muscle  + 4*T1_muscle -T0_muscle + input(i)*s_muscle * abs(pressure_xr_muscle).^2)./3;
%      T2_blood=T1(index_tissueblood:index_bloodmuscle,:);   
    T2_blood(1,:)=[];T2_muscle(1,:)=[];
    T2=[T2_multi; T2_blood; T2_muscle];
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
figure();
surf(z,xx,pressure_xr');
shading interp;
axis equal;
colormap(jet);
figure ()
surf(T2'-273);
shading interp;
axis equal;
colormap(jet);
hold on;
zv1=repmat(221,1,497);
zv2=repmat(241,1,497);
plot3(zv1,1:497,T2(221,:)-273,'r');
plot3(zv2,1:497,T2(241,:)-273,'r');
% scatter(z_index,round(gridNum_r/2,T2(z_index,round(gridNum_r/2)),'w');
figure();
plot(T_focus);
figure();
contourf(ThermalDose',[240,240]);

       