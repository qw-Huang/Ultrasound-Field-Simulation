%ʹ��dS1����Դ���ַ�������2dƽ��xz�棬�ӿ�Բ�׵Ļ�������ǿ���㣬�ȽϿ��ײ��֡�ʣ�ಿ�ֺ���������������ѹ��ֵ����λ
clear all;
close all;
f0=1e6;%����Ƶ�ʺ�������
P=100;
medium = set_medium('lossless');%������ʣ����㣺ˮ�������Ըĳɶ��set_layered_medium
lambda = medium.soundspeed/f0;%����=c/f
k=2*pi/lambda;%����
R = 1 * 5 * 2 * lambda;%ROC���ʰ뾶
a = 1 * 5 * lambda;%ע�������a�ǿ׾���һ��
fnumber=R/(2*a);%����f-number=���ʰ뾶/�׾���2*a��
d = sqrt(R^2 - a^2);%���۽��㵽�׾����ĵľ���
u=normal_velocity(P,R,a,0,medium.density,medium.soundspeed);

%���������
xmin=-a;%�۲������ķ�Χ
xmax=-xmin;
ymax=0;
ymin=0;
zDiff=0.7*d;
zmin=R-zDiff;
zmax=R+zDiff;

dx = lambda/6; %�����Ĳ���
dz = lambda/6;

x=xmin:dx:xmax;%�����ķֲ�
z=zmin:dz:zmax;

% ��y=0ʱxzƽ��
y=0;

%�ӿ���ʣ�ಿ��
hole_a=15e-4;
pr_res=hole_rayleigh2(R,a,u,hole_a,x,z); %������ز����õ���ѹp

%���ڵ��Ŀ���
pr_hole=hole_rayleigh2(R,hole_a,u,0,x,z); %������ز����õ���ѹp

%����������
pr_all=hole_rayleigh2(R,a,u,0,x,z); %������ز����õ���ѹp

% %ѡ�����Ȥ�������
% [index]=find(I_pr>=0.25*I_pr_max);
% contribution_res=I_pr_res./I_pr;
% contribution_hole=I_pr_hole./I_pr;
% res_zeros_xz(index)=contribution_res(index);
% hole_zeros_xz(index)=contribution_hole(index);

%�Ƚ���ѹ�ķ�ֵ����λ
pr_all_amplitude=abs(pr_all);
pr_all_phase=angle(pr_all);
pr_sum=pr_hole+pr_res;
pr_sum_amplitude=abs(pr_sum);
pr_sum_phase=angle(pr_sum);
error_phase=abs(pr_sum_phase-pr_all_phase)./pr_all_phase;
error_amplitude=abs(pr_sum_amplitude-pr_all_amplitude)./pr_all_amplitude;
I_all=abs(pr_all).^2./(2*medium.soundspeed*medium.density);
I_hole=abs(pr_hole).^2./(2*medium.soundspeed*medium.density);
I_res=abs(pr_res).^2./(2*medium.soundspeed*medium.density);
max_hole=find_maxpoint(I_hole);
max_res=find_maxpoint(I_res);
max_all=find_maxpoint(I_all);
% %��ͼ
% % �Ա�ʣ�ಿ�ֺ�ȥ�����ֶ������ֲ��Ĺ���
figure(1);
surf(z*1000, x*1000,pr_all_amplitude); 
axis equal;
shading interp;
colorbar;
xlabel('z (mm) ');
ylabel('x (mm) ');
title('pr_all_amplitude(R=15mm,a=7.5mm,20%a)');
figure(2);
surf(z*1000, x*1000,abs(pr_hole));
axis equal;
shading interp;
colorbar;
xlabel('z (mm) ');
ylabel('x (mm) ');
title('pr_hole_amplitude(R=15mm,a=7.5mm,20%a)');
figure(3);
surf(z*1000, x*1000,abs(pr_res));
axis equal;
shading interp;
colorbar;xlabel('z (mm) ');
ylabel('x (mm) ');
title('pr_res_amplitude(R=15mm,a=7.5mm,20%a)');
 figure(4);
surf(z*1000, x*1000,error_amplitude); 
axis equal;
shading interp;
colorbar;
xlabel('z (mm) ');
ylabel('x (mm) ');
title('the error of amplitude(R=15mm,a=7.5mm,20%a)');
figure(5);
surf(z*1000, x*1000,error_phase);
axis equal;
shading interp;
colorbar;
xlabel('z (mm) ');
ylabel('x (mm) ');
title('the error of phase(R=15mm,a=7.5mm,20%a)');

% figure(6);
% plot(abs(pr_all(31,:)));
% hold on;
% plot(abs(pr_res(31,:)));
% hold on;
% plot(abs(pr_hole(31,:)));
% figure(7);
% plot(angle(pr_all(31,:)));
% hold on;
% plot(angle(pr_res(31,:)));
% hold on;
% plot(angle(pr_hole(31,:)));

% shading interp ;%ȥ������ƽ������
% axis equal;
% colorbar  %����ɫ��
% xlabel('z (mm) ');
% ylabel('x (mm) ');
% title('(R=15mm,a=7.5mm,30%a)');
% save('R15_hole0.3a_holeVSresidue.mat','R','a','pr','pr_hole','pr_res','I_pr_hole','I_pr_res','I_pr','res_zeros_xz','hole_zeros_xz');
