tic
% p_max=max(abs(p(:)));
% [index]=find_maxpoint(abs(p));
% xindex=index(1);
% yindex=index(2);
% zindex=index(3);
%  [focalsize_x,focalsize_y,focalsize_z] = focalsize(P3D,p_max,xi,yi,zi);%计算焦点大小
% f0=1e6;  % muscle;  % blood;
clear p_asa1  p_asa2 p_asa3 p_asa4 p_asa5 p0 p0_all p0_res p0_hole p_abs;

XStep=dx;
deltaz=dz;
delta_t = 0.5;
% P3D=p;%(:,yindex,:);

% Convert attenuation to Np/m
dBperNeper = 20 * log10(exp(1));
alphaskin=medium2.attenuationdBcmMHz/dBperNeper*100*(f0/1e6)^1.1;
alphafat=medium3.attenuationdBcmMHz/dBperNeper*100*(f0/1e6)^1.1;
alphamuscle=medium4.attenuationdBcmMHz/dBperNeper*100*(f0/1e6)^1.1;
alphafibroid=medium5.attenuationdBcmMHz/dBperNeper*100*(f0/1e6)^1.1;
Alpha_skin=alphaskin*ones(nx_new,ny_new,nz2);
Alpha_fat=alphafat*ones(nx_new,ny_new,nz3-1);
Alpha_muscle=alphamuscle*ones(nx_new,ny_new,nz4-1);
Alpha_fibroid=alphafibroid*ones(nx_new,ny_new,nz5-1);
Alpha=cat(3,Alpha_skin, Alpha_fat, Alpha_muscle, Alpha_fibroid);    %构造参数矩阵，不同组织参数不同
clear Alpha_skin Alpha_fat Alpha_muscle Alpha_fibroid;
rho_skin=medium2.density*ones(nx_new,ny_new,nz2);
rho_fat=medium3.density*ones(nx_new,ny_new,nz3-1);
rho_muscle=medium4.density*ones(nx_new,ny_new,nz4-1);
rho_fibroid=medium5.density*ones(nx_new,ny_new,nz5-1);
rho=cat(3,rho_skin, rho_fat, rho_muscle, rho_fibroid);
clear rho_skin rho_fat rho_muscle rho_fibroid;

c_skin=medium2.soundspeed*ones(nx_new,ny_new,nz2);
c_fat=medium3.soundspeed*ones(nx_new,ny_new,nz3-1);
c_muscle=medium4.soundspeed*ones(nx_new,ny_new,nz4-1);
c_fibroid=medium5.soundspeed*ones(nx_new,ny_new,nz5-1);
c=cat(3,c_skin, c_fat, c_muscle, c_fibroid);
clear c_skin c_fat c_muscle c_fibroid;

Q3D = Alpha.*abs(P3D).^2./(rho.*c);

kt_skin=medium2.thermalconductivity*ones(nx_new,ny_new,nz2);
kt_fat=medium3.thermalconductivity*ones(nx_new,ny_new,nz3-1);
kt_muscle=medium4.thermalconductivity*ones(nx_new,ny_new,nz4-1);
kt_fibroid=medium5.thermalconductivity*ones(nx_new,ny_new,nz5-1);
kt3D=cat(3,kt_skin, kt_fat, kt_muscle, kt_fibroid);           % kt =0.5;
kt = reshape(kt3D,size(kt3D,1)*size(kt3D,2)*size(kt3D,3),1);
clear kt_skin kt_fat kt_muscle kt_fibroid kt3D;

C_skin=medium2.specificheat*ones(nx_new,ny_new,nz2);
C_fat=medium3.specificheat*ones(nx_new,ny_new,nz3-1);
C_muscle=medium4.specificheat*ones(nx_new,ny_new,nz4-1);
C_fibroid=medium5.specificheat*ones(nx_new,ny_new,nz5-1);
C3D=cat(3,C_skin, C_fat, C_muscle, C_fibroid);              % C =3800;
C = reshape(C3D,size(C3D,1)*size(C3D,2)*size(C3D,3),1);
clear C_skin C_fat C_muscle C_fibroid C3D;

rhore = reshape(rho,size(rho,1)*size(rho,2)*size(rho,3),1);

a1 = rhore.*C/delta_t + 2*kt/(XStep*XStep)+kt/(deltaz*deltaz);
b1 = -kt/(2*XStep*XStep);
c1 = -kt/(2*deltaz*deltaz);
a2 = rhore.*C/delta_t - 2*kt/(XStep*XStep)-kt/(deltaz*deltaz);
b2 = -b1;
c2 = -c1;
 e = ones(size(P3D,1)*size(P3D,2)*size(P3D,3),1);

 Qc =0;% Wb * rho_b * Cb;
nx=nx_new;
ny=ny_new;
A = spdiags([c1,b1,b1,a1,b1,b1,c1],[-(nx*ny),-nx,-1,0,1,nx,(nx*ny)],length(e),length(e));

B = spdiags([c2,b2,b2,a2,b2,b2,c2],[-(nx*ny),-nx,-1,0,1,nx,(nx*ny)],length(e),length(e));

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

% figure(1);
% surf(squeeze(Trise(:,81,50:177)));
% % surf(squeeze(Trise));
% shading interp;
% colormap(jet);
% toc



% error=abs(Tnew-Trise_xy)./Tnew;
% surf(error);
% shading interp;
% axis equal;
% colormap(jet);
