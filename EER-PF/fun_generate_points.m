function [x,step,angel] = fun_generate_points(num,length,width,step_size) 

x=zeros(num,2);         %100个位置
step=zeros(num-1,1);    %每一次步长
angel=zeros(num-1,1);   %每一次角度
% 
x(1,1)=unifrnd(0,width);%初始位置
x(1,2)=unifrnd(0,length);

%步长和角度
for i = 1:num-1
    angel(i) = rand(1,1)*360;
    step(i) = step_size;
end

%轨迹
for i=2:num  %真实位置信息
    x(i,1)= x(i-1,1)+step(i-1)*cos(angel(i-1)/180*pi);
    x(i,2)= x(i-1,2)+step(i-1)*sin(angel(i-1)/180*pi);
    %确保在场景内
    while x(i,1)>width || x(i,1)<0||x(i,2)>length || x(i,2)<0%如果出界，就选择向反方向走
    	angel(i-1)=angel(i-1) + 90;
        x(i,1)= x(i-1,1)+step(i-1)*cos(angel(i-1)/180*pi);
        x(i,2)= x(i-1,2)+step(i-1)*sin(angel(i-1)/180*pi);
    end
end
