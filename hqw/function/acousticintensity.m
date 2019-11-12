function [I] = acousticintensity(pressure,density,soundspeed)
%输入参数为声压、密度、声速，声强计算 I=|p|的平方/(ρc)，输出声强的值
I=abs(pressure).*abs(pressure)/(density*soundspeed); 
end

