
%�ı�f-number��������ѧ���Եı仯
clc;
clear all;
close all;
tic 

%�������е���ѹ��ֵ��ʽ�����󵼣��õ����㴦����һ�׵�������0��ʱ��õ�������ѹ���ֵ��Ӧ��w

n=0:0.005:0.5;

f0=1e6;%����Ƶ�ʺ�������
medium = set_medium('lossless');%������ʣ����㣺ˮ��
lambda = medium.soundspeed/f0;%����=c/f
P=100;
A=75e-3;
k=2*pi/lambda;%����
% focus0=[4e-3 14.9e-3 14.9e-3];
syms w;
a=30e-3;
%���������й�ʽ������ǰ�ƾ����濪�ױ����ı仯
for i=1:length(n)
    hole_a=n(i)*a;
     u(i)=normal_velocity(P,A,a,hole_a,medium.density,medium.soundspeed);
    y2=asin(a/A);
    y1=asin(hole_a/A);
    U2=sqrt(w.^2+2*A.*(A-w).*(1-cos(y2)));
    U1=sqrt(w.^2+2*A.*(A-w).*(1-cos(y1)));
    p_abs=2*medium.density*medium.soundspeed*u(i)*A./abs(A-w).*abs(sin(k.*(U2-U1)./2));
    dp=diff(p_abs,w); %�Ժ���p(w)�󵼣���ROC����һ�׵���Ϊ0�ĵ������ֵ���w����
    w0=vpasolve(dp==0,w,74.9e-3);
    w_forward(i)=A-double(w0); %��ROC-w���õ�ǰ�ƾ���    
end
toc
scatter(n,w_forward*1000,5,'.','b');
hold on;

a=37.5e-3;
%���������й�ʽ������ǰ�ƾ����濪�ױ����ı仯
for i=1:length(n)
    hole_a=n(i)*a;
     u(i)=normal_velocity(P,A,a,hole_a,medium.density,medium.soundspeed);
    y2=asin(a/A);
    y1=asin(hole_a/A);
    U2=sqrt(w.^2+2*A.*(A-w).*(1-cos(y2)));
    U1=sqrt(w.^2+2*A.*(A-w).*(1-cos(y1)));
    p_abs=2*medium.density*medium.soundspeed*u(i)*A./abs(A-w).*abs(sin(k.*(U2-U1)./2));
    dp=diff(p_abs,w); %�Ժ���p(w)�󵼣���ROC����һ�׵���Ϊ0�ĵ������ֵ���w����
    w0=vpasolve(dp==0,w,74.9e-3);
    w_forward(i)=A-double(w0); %��ROC-w���õ�ǰ�ƾ���    
end
toc
scatter(n,w_forward*1000,5,'.','r');
hold on;

a=45e-3;
%���������й�ʽ������ǰ�ƾ����濪�ױ����ı仯
for i=1:length(n)
    hole_a=n(i)*a;
     u(i)=normal_velocity(P,A,a,hole_a,medium.density,medium.soundspeed);
    y2=asin(a/A);
    y1=asin(hole_a/A);
    U2=sqrt(w.^2+2*A.*(A-w).*(1-cos(y2)));
    U1=sqrt(w.^2+2*A.*(A-w).*(1-cos(y1)));
    p_abs=2*medium.density*medium.soundspeed*u(i)*A./abs(A-w).*abs(sin(k.*(U2-U1)./2));
    dp=diff(p_abs,w); %�Ժ���p(w)�󵼣���ROC����һ�׵���Ϊ0�ĵ������ֵ���w����
    w0=vpasolve(dp==0,w,74.9e-3);
    w_forward(i)=A-double(w0); %��ROC-w���õ�ǰ�ƾ���    
end
toc
scatter(n,w_forward*1000,5,'.','k');
xlabel('hole(a)/a');
ylabel('focus forward distance(mm)');

% % % %���ݹ�ʽ�� ����-6dB��Χ�仯��Χ
% 
% for i=1:length(n)
%     w=0:0.000001:22e-3;
%     hole_a=n(i)*a(ia);
%     y2=asin(a(ia)/A);
%     y1=asin(hole_a/A);
%     U2=sqrt(w.^2+2*A.*(A-w).*(1-cos(y2)));
%     U1=sqrt(w.^2+2*A.*(A-w).*(1-cos(y1)));
%     p_abs=2*medium.density*medium.soundspeed*u*A./abs(A-w).*abs(sin(k.*(U2-U1)./2));
%     p_max(ia,i)=max(p_abs(:));
%     index_z=find(p_abs>=max(p_max(ia,i)*0.5));
%     index_z(ia,:)=index_z;
%     dB(ia,i)=w(max(index_z))-w(min(index_z));
% end

% figure(2);
% scatter(n,dB(1,:),'.','b');
% xlabel('���״�С');
% ylabel('����-6dB��Χ�仯');
% hold on;
% scatter(n,dB(2,:),'.','r');
% scatter(n,dB(3,:),'.','k');
% 
% figure(3);
% scatter(n,p_max(1,:),'.','b');
% xlabel('���״�С');
% ylabel('pmax');
% hold on;
% scatter(n,p_max(2,:),'.','r');
% scatter(n,p_max(3,:),'.','k');
