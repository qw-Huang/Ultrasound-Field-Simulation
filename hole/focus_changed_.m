%ʹ���Զ���rayleigh���ּ��㿪�׶Խ�����ǿ��-6dB������ǰ�Ƶ�Ӱ�� ��ֵ��� ������� ���ں����۽�����бȽ�
clc;
clear all;
ROC=15e-3;
a=7.5e-3;
u=1;
n=0:0.001:0.5;
x=0;
z=8e-3:0.000001:22e-3;
%����ǰ�ƺͿ��״�С�Ĺ�ϵ
for i=1:length(n)
   hole_a(i)=n(i)*a;
   pr_axial=hole_rayleigh2(ROC,a,u,hole_a(i),x,z);
   pr_axial=abs(pr_axial);
   pr_max(i)=max(pr_axial(:));
   index_focus_z(i)=find(pr_axial==pr_max(i));
   focus_forward(i)=ROC-z(index_focus_z(i));
    
   %����-6dB�Ϳ��״�С�Ĺ�ϵ
   index_dB_axial=find(pr_axial>=0.5*pr_max(i));
   index_axial_forward(i)=min(index_dB_axial);
   index_axial_backward(i)=max(index_dB_axial);
   axial_dB(i)=z(index_axial_backward(i))-z(index_axial_forward(i));
end
figure(1);
scatter(n,focus_forward,'.','k');
figure(2);
scatter(n,pr_max,'.','k');
figure(3);
scatter(n,axial_dB,'.','k');


% %����-6dB�Ϳ��״�С�Ĺ�ϵ
% x=-7.5e-3::0:
% for i=1:length(n)
%    hole_a(i)=n(i)*a;
%    pr_radial=hole_rayleigh2(ROC,a,u,hole_a(i),x,z(index_focus_z(i)));
%    pr_radial=abs(pr_radial);
%    index_x(i)=find(pr_radial>=0.5*pr_max);
%    focus_forward(i)=ROC-z(index_focus_z(i));
% end