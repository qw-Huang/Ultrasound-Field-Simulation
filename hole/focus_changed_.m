%使用自定义rayleigh积分计算开孔对焦点声强、-6dB、焦点前移的影响 数值结果 进行拟合 便于和理论结果进行比较
clc;
clear all;
ROC=15e-3;
a=7.5e-3;
u=1;
n=0:0.001:0.5;
x=0;
z=8e-3:0.000001:22e-3;
%焦点前移和开孔大小的关系
for i=1:length(n)
   hole_a(i)=n(i)*a;
   pr_axial=hole_rayleigh2(ROC,a,u,hole_a(i),x,z);
   pr_axial=abs(pr_axial);
   pr_max(i)=max(pr_axial(:));
   index_focus_z(i)=find(pr_axial==pr_max(i));
   focus_forward(i)=ROC-z(index_focus_z(i));
    
   %轴向-6dB和开孔大小的关系
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


% %径向-6dB和开孔大小的关系
% x=-7.5e-3::0:
% for i=1:length(n)
%    hole_a(i)=n(i)*a;
%    pr_radial=hole_rayleigh2(ROC,a,u,hole_a(i),x,z(index_focus_z(i)));
%    pr_radial=abs(pr_radial);
%    index_x(i)=find(pr_radial>=0.5*pr_max);
%    focus_forward(i)=ROC-z(index_focus_z(i));
% end