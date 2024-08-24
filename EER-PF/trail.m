 %由randomwalk生成轨迹
  function [RMSE1,RMSE2,RMSE3,RMSE4,variance1,variance2,variance3,variance4] = trail(~)
%% 全局变量
global x xpr xarr angel
step_mean = 0.1;
step_sigma = 0.1;
angel_mean = 10;%5
angel_sigma = 0.1;%0.25
sum1 = 0;
sum2 = 0;
sum3 = 0;
sum4 = 0;
N = 6000;%粒子数量
European_distance1 = zeros(10,1);
European_distance2 = zeros(10,1);
European_distance3 = zeros(10,1);
European_distance4 = zeros(10,1);
% Cosine_distance1 = zeros(10,1);
% Cosine_distance2 = zeros(10,1);
% %%
num       =100;	 %行走次数
length    =50;  %场景长
width     =50;  %场景宽
step_size =2;   %根据场景决定吧
%%
%生成实际轨迹与画图
[x,step,angel] =fun_generate_points(num,length,width,step_size);
%%
%生成预测轨迹与画图
[xmen,xpr,sr_xpr,syr_xpr,xarr] = fun_pr_trail(x,num,step,angel,step_mean,step_sigma,angel_mean,angel_sigma,N);
% 粒子均值
%% 欧式距离
for i = 1:num
    European_distance1(i) = sqrt((xmen(i,1) - x(i,1))^2 + (xmen(i,2) - x(i,2))^2);
    European_distance2(i) = sqrt((xpr(i,1) - x(i,1))^2 + (xpr(i,2) - x(i,2))^2);
    European_distance3(i) = sqrt((sr_xpr(i,1) - x(i,1))^2 + (sr_xpr(i,2) - x(i,2))^2);
    European_distance4(i) = sqrt((syr_xpr(i,1) - x(i,1))^2 + (syr_xpr(i,2) - x(i,2))^2);
    if European_distance1(i) < European_distance2(i)
        disp(i);
    end
end
distance1 = mean(European_distance1);
distance2 = mean(European_distance2);
distance3 = mean(European_distance3);
distance4 = mean(European_distance4);

variance1 = var(European_distance1);
variance2 = var(European_distance2);
variance3 = var(European_distance3);
variance4 = var(European_distance4);

for i = 1:num
    sum1 = sum1 + (European_distance1(i) - distance1)^2;
    sum2 = sum2 + (European_distance2(i) - distance2)^2;
    sum3 = sum3 + (European_distance3(i) - distance3)^2;
    sum4 = sum4 + (European_distance4(i) - distance4)^2;
end
RMSE1 = sqrt(sum1/num);
RMSE2 = sqrt(sum2/num);
RMSE3 = sqrt(sum3/num);
RMSE4 = sqrt(sum4/num);