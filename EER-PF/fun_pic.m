function [] = fun_pic(x,num,length,width,linecolor) 

pointcolor='r.';
ifdrawpoint=0;
%开始
plot(x(1,1),x(1,2),'rp','markersize',10,'linewidth',2)
text(x(1,1),x(1,2),'begin');
hold on;
%结束
plot(x(num,1),x(num,2),'rp','markersize',10,'linewidth',3)
plot([x(num-1,1) x(num,1)],[x(num-1,2) x(num,2)],linecolor,'linewidth',2);
text(x(num,1),x(num,2),'end');

%之间
for i=2:num-1
    if ifdrawpoint==1
        plot(x(i,1),x(i,2),pointcolor,'markersize',10,'linewidth',2)
    end
    plot([x(i-1,1) x(i,1)],[x(i-1,2) x(i,2)],linecolor,'linewidth',2);%只有线段
end

%边界，因为imu轨迹会出框，所以画一下
rectangle('position',[0 0 width length] );

axis equal
axis([-5 width+5 -5 length+5])
hold off;
%想办法画出箭头