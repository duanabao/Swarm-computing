function [x,step,angel] = fun_generate_points(num,length,width,step_size) 

x=zeros(num,2);         %100��λ��
step=zeros(num-1,1);    %ÿһ�β���
angel=zeros(num-1,1);   %ÿһ�νǶ�
% 
x(1,1)=unifrnd(0,width);%��ʼλ��
x(1,2)=unifrnd(0,length);

%�����ͽǶ�
for i = 1:num-1
    angel(i) = rand(1,1)*360;
    step(i) = step_size;
end

%�켣
for i=2:num  %��ʵλ����Ϣ
    x(i,1)= x(i-1,1)+step(i-1)*cos(angel(i-1)/180*pi);
    x(i,2)= x(i-1,2)+step(i-1)*sin(angel(i-1)/180*pi);
    %ȷ���ڳ�����
    while x(i,1)>width || x(i,1)<0||x(i,2)>length || x(i,2)<0%������磬��ѡ���򷴷�����
    	angel(i-1)=angel(i-1) + 90;
        x(i,1)= x(i-1,1)+step(i-1)*cos(angel(i-1)/180*pi);
        x(i,2)= x(i-1,2)+step(i-1)*sin(angel(i-1)/180*pi);
    end
end
