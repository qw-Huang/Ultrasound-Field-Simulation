function [I] = acousticintensity(pressure,density,soundspeed)
%�������Ϊ��ѹ���ܶȡ����٣���ǿ���� I=|p|��ƽ��/(��c)�������ǿ��ֵ
I=abs(pressure).*abs(pressure)/(density*soundspeed); 
end

