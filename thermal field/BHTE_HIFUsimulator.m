% clear all
close all;
load('pressure_xz_muscle1MPa.mat');
rho1=1140;
c1=1570;
rho2=1140;
c2=1570;
alpha=9;
K=1;%谐波次数
r=x;%x_new;%((1+nx)/2:nx);
% z=z_new;
J=length(r);
M=length(z);
P3D=p;
pressure=P3D(:,127,:);%((1+nx)/2:nx,:);
H = alpha*abs(pressure).^2/rho1/c1;  %/d;
% material 1:
C1 = 3480;	% heat capacity (J/kg/K)
k1 = 0.4980;	% thermal conductivity (W/m/K)
w1 = 0;		% perfusion rate (kg/m^3/K)
C2 = 3480;	% heat capacity (J/kg/K)
k2 = 0.4980;	% thermal conductivity (W/m/K)
w2 = 0;		% perfusion rate (kg/m^3/K)

% ambient temperature:
T0 = 37;	% (degrees C)

dt=0.5;     %(unit:s)
% N = 61;		% total number of integration steps
% Ts= dt*(N-1);			% total simulation duration
% t=0:dt:Ts;
% input=1;
D1 = k1/C1/rho1;			% diffusivity of material 1
D2 = k2/C2/rho2;			% diffusivity of material 2
P1 = w1/rho1;
P2 = w2/rho2;
% z_/d/Z=1 介质1传递距离的网格点数

% continuous sonication
t_i = 50.5;	% initial sonication duration (s)

% % periodic sonications
n_c = 0;		% number of pulse cycles
% D = 20;		% duty cycle (%)
% t_p = 0.5;	% pulse cycle period (s)
% cooling period
t_c = 0;	% cool-off duration (s)

% determine number of integration steps at each stage:
N_i = round(t_i/dt);
if(n_c==0)
  N_p=0;
else
  N_p = round(t_p/dt);
end
N_c = round(t_c/dt);
N = N_i+n_c*N_p+N_c;		% total number of integration steps
Ts = dt*(N-1);			% total simulation duration
t=0:dt:Ts;
% build input vector, which catalogs the HIFU beam on/off operation 
if(N_p~=0)
    pulse = zeros(1,N_p);
  for n=1:round(0.01*D*N_p)
    pulse(n) = 1;
  end
  pulses = pulse;
  for m=2:n_c
    pulses = [pulses,pulse];
  end
end
if(N_c~=0)
  cooloff = zeros(1,N_c);
end 
if(N_i~=0)
  initial = ones(1,N_i);
end
for n=1:N_i
  input(n) = initial(n);
end
for n=N_i+1:N_i+n_c*N_p
  input(n) = pulses(n-N_i);
end
for n=N_i+n_c*N_p+1:N
    input(n) = cooloff(n-N_i-n_c*N_p);
end

% create Crank-Nicolson operators for BHT:
[A,B] = BHT_operators(r,z,D1,D2,P1,P2,dt,1);
m_ = round(M);			% index @ material interface 
if(m_>M) m_ = M; end
H(:,1:m_) = H(:,1:m_)/C1/rho1;	% rescale H to degrees/second
H(:,m_+1:M) = H(:,m_+1:M)/C2/rho2;

JM = J*M;
Q = zeros(JM,1);			% heating source vector
T = zeros(JM,1);			% temperature vector
D = zeros(JM,1);			% thermal dose vector
Tmat = zeros(J,M);			% temp matrix
Dmat = zeros(J,M);			% dose matrix
Q = vektorize(Q,H,J,M);			% column-stack the heat matrix
s1 = zeros(JM,1);			% slopes for IRK integrator
s2 = zeros(JM,1);
Tpeak = zeros(N,1);			% temp vs time vector
Tmax_vec = zeros(JM,1);			% max temp distribution vector
Tmax_mat = zeros(J,M);			% max temp distribution matrix
tt = 0;

% Integrate BHT:
fprintf('\tIntegrating BHT equation...\n')
fprintf('\t\tt (sec)\ttime (hr:min:sec)\n')
t_start = clock;
p1 = 0;
for n=1:N-1
  s1 = A*T + input(n)*Q;
  s2 = B \ (A*(T+0.5*dt*s1) + input(n+1)*Q);
  T = T + 0.5*dt*(s1+s2);
%   D = D + equivalent_time(T,JM,T0);
   Tpeak(n+1) = max(T);
   T_focus(n)=Tpeak(n+1); 
%   if(Tpeak(n+1)>tt)
%     tt = Tpeak(n+1);
    Tmax_vec = T;
%     t_peak = t(n+1);
%   end
  p2 = floor(10*(n+1)/N);
  p1 = timing(p1,p2,t_start,t(n+1),1);
end

Tmax_mat = matrixize(Tmax_vec,Tmax_mat,J,M);
Tnew=Tmax_mat+37;
figure();
surf(Tmax_mat);
shading interp;
axis equal;
colormap(jet);
figure();
plot(t(2:length(t)),T_focus+37);
% 
load('Thermaldose_muscle1Mpa.mat')
% Trise=squeeze(Trise(:,164,:)+37);
error=abs(Tnew-squeeze(Trise))./Trise;
figure();
surf(error);
shading interp;
axis equal;
colormap(jet);
max(error(:))
