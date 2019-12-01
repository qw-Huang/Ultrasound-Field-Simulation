function [output] = hole_ring_r(R,a,hole_ring_back,hole_ring_after,x,z)
%HOLE_RING �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%ʹ��dS1����Դ���ַ�������2dƽ��xz�棬�ӿ��׵Ļ�������ǿ����

f0=1e6;%����Ƶ�ʺ�������
P=100;
medium = set_medium('lossless');%������ʣ����㣺ˮ�������Ըĳɶ��set_layered_medium
lambda = medium.soundspeed/f0;%����=c/f
k=2*pi/lambda;%����
fnumber=R/(2*a);%����f-number=���ʰ뾶/�׾���2*a��
d = sqrt(R^2 - a^2);%���۽��㵽�׾����ĵľ���
u=normal_velocity(P,R,a,0,medium.density,medium.soundspeed);

nx=length(x);%�����ĵ���
nz=length(z);

error_zeros_xz=zeros(nx,nz); %����ȫ��������

%��������ɢ��Ϊ����Դ  ע����ɢ��֮�����Դ�Ĵ�С
dr=lambda/6;%����ָ�ɺܶ��Բ������Ӧ�İ뾶��r��dr�ǰ뾶���ӵĲ���
ndr=round((a-0)/dr);
dr=(a-0)/ndr;
r_back=0:dr:a-dr;%����Դǰһ�λ�����Ӧ��r
r_after=dr:dr:a;%����Դ��һ�λ�����Ӧ��r
r=r_after-dr/2;%��i��������Ӧ�����ĵ��r
Sm=lambda/6;%����㻷����2pi*a����ɢ�����Ӧ��dS�Ļ���С��lambda/6
median=fix(length(r)/2)+1;%fix��λȡ��
ntheta=round(2*pi*r(median)/Sm);%�����һ�������Ļ��ֵ���,ȡ��
dtheta=2*pi./ntheta;%����ȡ�����µ���ÿ������Դ�Ķ�Ӧ����
theta_after=dtheta:dtheta:2*pi;%ÿ��������ɢ�ɶ�����Ӧ�Ļ�������
theta_back=0:dtheta:(2*pi-dtheta);%ÿ��������ɢ�ɶ�����Ӧ�Ļ�������
theta=theta_after-dtheta/2;%��i��������ɢ�ɵ���Դ��dS�м�ĵ��Ӧ�Ļ�������
X=cos(theta)'*r;%��Դ��x����  ������ʽ ntheta*nr
Y=sin(theta)'*r;%��Դ��y����
Z0=R-sqrt(R*R-r.*r);%����Դ��ά�ռ��е�z����
Z=repmat(Z0,ntheta,1); % repmat( A , m , n )���������������ڴ�ֱ������m�Σ���ˮƽ������n�Ρ�

%���廻������������
u_ring=repmat(u,ntheta,ndr);
back_index=find(r_back<=hole_ring_back,1,'last');%�жϻ����ھ������ĸ���Χ��
after_index=find(r_after>=hole_ring_after,1,'first');%�жϻ����⾶�����ĸ���Χ��
if back_index==after_index
    index_ring=back_index;
else index_ring=back_index:after_index;
end
u_ring(:,index_ring)=[0]; %�����ײ��ַ���������Ϊ0

% ��y=0ʱxzƽ��
y=0;

%�󿪿�Ϊ���ε���ѹ
tic
for ix=1:nx  %�۲������x���������
    for iz=1:nz  %�۲������z���������
        rn=sqrt((X-x(ix)).^2+(Y-y).^2+(Z-z(iz)).^2);%�۲�㵽��Դ�ľ���
        dS_ring=r.*dtheta.*R.*(asin(r_after./R)-asin(r_back./R));%ÿһ������ɢ��dS�Ĵ�С����i����ɢ�ĵ���Դ���dS=��i������ɢ����*dr��Ӧ�Ķ̻�
        dS=repmat(dS_ring,ntheta,1); % repmat( A , m , n )���������������ڴ�ֱ������m�Σ���ˮƽ������n�Ρ�
        A=dS.*exp(-1i.*k.*rn)./rn.*u_ring;
        B=sum(sum(A));%��������õ�ֵ�ۼ�
        pr_hole_ring_res(ix,iz)=1i*medium.density*medium.soundspeed*k/(2*pi)*B; %������ز����õ���ѹp
    end
end
toc
output = pr_hole_ring_res;
end

