
%�Ա������Ƶ��������ɢ��������������ѧ���Եı仯
clc;
clear all;
close all;
tic 

%�������е���ѹ��ֵ��ʽ�����󵼣��õ����㴦����һ�׵�������0��ʱ��õ�������ѹ���ֵ��Ӧ��w
syms w;
n=0:0.002:0.5;
a=7.5e-3;
f0=1e6;%����Ƶ�ʺ�������
medium = set_medium('lossless');%������ʣ����㣺ˮ��
lambda = medium.soundspeed/f0;%����=c/f
u=1;
A=15e-3;
ROC=15e-3;
k=2*pi/lambda;%����

%���������й�ʽ������ǰ�ƾ����濪�ױ����ı仯
for i=1:length(n)
    hole_a=n(i)*a;
    y2=asin(a/A);
    y1=asin(hole_a/A);
    U2=sqrt(w.^2+2*A.*(A-w).*(1-cos(y2)));
    U1=sqrt(w.^2+2*A.*(A-w).*(1-cos(y1)));
    p_abs=2*medium.density*medium.soundspeed*u*A./abs(A-w).*abs(sin(k.*(U2-U1)./2));
    dp=diff(p_abs,w); %�Ժ���p(w)�󵼣���ROC����һ�׵���Ϊ0�ĵ������ֵ���w����
    w0=vpasolve(dp==0,w,14.9e-3);
    w_forward(i)=A-double(w0); %��ROC-w���õ�ǰ�ƾ���    
end
% figure(1);
% scatter(n,w_forward,'.','k');
% xlabel('���״�С');
% ylabel('ǰ�ƾ���');
% hold on;
toc

x=0;
z=0:0.000001:22e-3;
%����rayleigh���ּ���ǰ�ơ�������ѹ������-6dB
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
% figure(1);
% scatter(n,focus_forward,'.','b');
% figure(2);
% scatter(n,pr_max,'.','b');
% hold on;
% figure(3);
% scatter(n,axial_dB,'.','b');
% 

% % %���ݹ�ʽ�� ����-6dB��Χ�仯��Χ
w=0:0.000001:22e-3;
for i=1:length(n)
    hole_a=n(i)*a;
    y2=asin(a/A);
    y1=asin(hole_a/A);
    U2=sqrt(w.^2+2*A.*(A-w).*(1-cos(y2)));
    U1=sqrt(w.^2+2*A.*(A-w).*(1-cos(y1)));
    p_abs=2*medium.density*medium.soundspeed*u*A./abs(A-w).*abs(sin(k.*(U2-U1)./2));
    p_max(i)=max(p_abs(:));
    index_z=find(p_abs>=max(p_max(i)*0.5));
    dB(i)=w(max(index_z))-w(min(index_z));
end
% figure(3);
% scatter(n,dB,'.','k');
% xlabel('���״�С');
% ylabel('����-6dB��Χ�仯');
% figure(2);
% scatter(n,p_max,'.','k');
% xlabel('���״�С');
% ylabel('pmax');

%��ͼ
figure(1);
scatter(n,w_forward*1000,5,'.','k');
xlabel('hole(a)/a');
ylabel('focus forward distance(mm)');
hold on;
scatter(n,focus_forward*1000,5,'.','b');

figure(2);
scatter(n,dB*1000,5,'.','k');
xlabel('hole(a)/a');
ylabel('axial distance -6dB(mm)');
hold on;
scatter(n,axial_dB*1000,5,'.','b');

figure(3);
scatter(n,p_max/10^6,5,'.','k');
xlabel('hole(a)/a');
ylabel('pmax(MPa)');
hold on;
scatter(n,pr_max/10^6,5,'.','b');


% 
% f_number=0.51:0.05:5;
% hole_a=0;
% %�ı�f-number ����������ǰ�����
% for i=1:length(f_number)
%     a=A/(2*f_number(i));
%     y2=asin(a/A);
%     y1=asin(hole_a/A);
%     U2=sqrt(w.^2+2*A.*(A-w).*(1-cos(y2)));
%     U1=sqrt(w.^2+2*A.*(A-w).*(1-cos(y1)));
%     p=2*medium.density*medium.soundspeed*u*A./abs(A-w).*abs(sin(k.*(U2-U1)./2));
%     index=find(p==max(p(:)));
%     focus_forward(i)=A-w(index);
% end
% figure(3);
% scatter(f_number,focus_forward,'.','k');
% xlabel('f-number');
% ylabel('ǰ�ƾ���');
% toc