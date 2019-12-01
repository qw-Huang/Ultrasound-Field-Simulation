%����ʱ�� ����ά�Ĵ�����иĽ���ֻ�����ķ�֮һ������ 
close all;
clear all;
clear all;
profile on;
f0=1e6;%����Ƶ�ʺͷ�������
P=100;
medium = set_medium('lossless');%������ʣ�����->ˮ���ɸĳɶ�㣬�ú���set_layered_medium
lambda = medium.soundspeed/f0;%����=c/f
k=2*pi/lambda;%����

R =10* 5 * 2 * lambda;%ROC���ʰ뾶
a =10* 5 * lambda;%ע������a�ǿ׾���һ��
fnumber=R/(2*a);%����f-number=���ʰ뾶/�׾���2*a��
d = sqrt(R^2 - a^2);%���۽��㵽�׾����ĵľ���
u=normal_velocity(P,R,a,0,medium.density,medium.soundspeed);

%���������
xmin=-0.005;%�۲������ķ�Χ
xmax=0;
ymin=-0.005;
ymax=0;
zDiff=0.01;
zmin=R-zDiff;
zmax=R+zDiff;

dx = lambda/6; %�����Ĳ���
dy = lambda/6; 
dz = lambda/6;

x=xmin:dx:xmax;%�����ķֲ�
y=ymin:dy:ymax;
z=zmin:dz:zmax;

nx=length(x);%�����ĵ���
ny=length(y);
nz=length(z);

%��������ɢ��Ϊ����Դ  ע����ɢ��֮�����Դ�Ĵ�С
dr=lambda/6;%����ָ�ɺܶ��Բ������Ӧ�İ뾶��r��dr�ǰ뾶���ӵĲ���
r_back=0:dr:a-dr;%����Դǰһ�λ�����Ӧ��r
r_after=dr:dr:a;%����Դ��һ�λ�����Ӧ��r
r=r_after-dr/2;%��i��������Ӧ�����ĵ��r
Sm=lambda/6;%����㻷����2pi*a����ɢ�����Ӧ��dS�Ļ���С��lambda/6
ntheta=round(2*pi*(a-dr/2)/Sm);%�����һ�������Ļ��ֵ���,ȡ��
dtheta=2*pi./ntheta;%����ȡ�����µ���ÿ������Դ�Ķ�Ӧ����
theta_after=dtheta:dtheta:2*pi;%ÿ��������ɢ�ɶ�����Ӧ�Ļ�������
theta_back=0:dtheta:(2*pi-dtheta);%ÿ��������ɢ�ɶ�����Ӧ�Ļ�������
theta=theta_after-dtheta/2;%��i��������ɢ�ɵ���Դ��dS�м�ĵ��Ӧ�Ļ�������
X=cos(theta)'*r;%��Դ��x����  ������ʽ ntheta*nr
Y=sin(theta)'*r;%��Դ��y����
Z_ring=R-sqrt(R*R-r.*r);%����Դ��ά�ռ��е�z����
Z=repmat(Z_ring,ntheta,1); % repmat( A , m , n )���������������ڴ�ֱ������m�Σ���ˮƽ������n�Ρ�

%rayleigh���ּ�����ά�ռ������  
tic
for ix=1:nx  %�۲������x��������껮��
    for iy=1:ny %�۲������y��������껮��
       for iz=1:nz  %�۲������z��������껮��
            rn=sqrt((X-x(ix)).^2+(Y-y(iy)).^2+(Z-z(iz)).^2);%�۲�㵽��Դ�ľ���
            dS_ring=r.*dtheta.*R.*(asin(r_after./R)-asin(r_back./R));%ÿһ������ɢ��dS�Ĵ�С����i����ɢ�ĵ���Դ���dS=��i������ɢ����*dr��Ӧ�Ķ̻�
            dS=repmat(dS_ring,ntheta,1); % repmat( A , m , n )���������������ڴ�ֱ������m�Σ���ˮƽ������n�Ρ�
            A=dS.*exp(-1i.*k.*rn)./rn;
            B=sum(sum(A));%��������õ�ֵ�ۼ�
            pr_1(ix,iy,iz)=1i*medium.density*u*medium.soundspeed*k/(2*pi)*B; %������ز����õ���ѹp   
        end
    end
end
toc
pr_3=flipud(pr_1);%����ѹ������з�ת���൱�ڹ���x=0��Գ�
pr_13=[pr_1;pr_3];%��������������ƴ�������������м�һ�лḴ������
[row,column]=size(pr_13);
row_median=row/2;%ȡ�м���
pr_13(row_median,:,:)=[]; %ɾ���м��ظ���һ��
pr_24=fliplr(pr_13);%����ѹ����������ҷ�ת���൱�ڹ���y=0��Գ�
pr_all=[pr_13,pr_24];%��������������ƴ�������������м�һ�лḴ������
[row,column,c]=size(pr_all);
col_median=column/2;%ȡ�м���
pr_all(:,col_median,:)=[]; %ɾ���м��ظ���һ��

%��ѹת��Ϊ��������
I_pr=acousticintensity(pr_all,medium.density,medium.soundspeed); 
I_pr_nor=I_pr./max(I_pr(:));
max_index=find_maxpoint(I_pr);%��������λ������
y_index=max_index(2);
x=-0.005:dx:0.005;
surf(z*1000,x*1000,squeeze(I_pr(:,y_index,:)));
shading interp;
axis equal;
profile viewer



