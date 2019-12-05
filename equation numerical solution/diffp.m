%�������е���ѹ��ֵ��ʽ�����󵼣��õ����㴦����һ�׵�������0��ʱ��õ�������ѹ���ֵ��Ӧ��w
clear all;
syms w;
% w=0:2.5e-4:14e-3;
f0=1e6;%����Ƶ�ʺ�������
medium = set_medium('lossless');%������ʣ����㣺ˮ��
lambda = medium.soundspeed/f0;%����=c/f
u=1;
A=15e-3;
k=2*pi/lambda;%����
y2=pi/5;
% y1=0.35*y2;
y1=0;
n=sin(y1)/sin(y2);
U2=sqrt(w.^2+2*A.*(A-w).*(1-cos(y2)));
U1=sqrt(w.^2+2*A.*(A-w).*(1-cos(y1)));
p=2*medium.density*medium.soundspeed*u*A./abs(A-w).*abs(sin(k.*(U2-U1)./2));
dp=diff(p,w);
w0=vpasolve(dp==0,w,14e-3);
w1=subs(w0);
% w=0:2.5e-4:14e-3;

