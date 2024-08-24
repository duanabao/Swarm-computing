function [syr_xx,syr_xpart,syrw] = SystematicResampling(syr_xarr,syr_true,N)
%%
%变量定义
syr_xpart = zeros(N,2);%粒子坐标信息
syr_weight = zeros(N,1);%粒子权重
R =1;
k = 0;
%%
%算法实现
for i = 1:N
    d = sqrt((syr_xarr(i,1) - syr_true(1))^2+(syr_xarr(i,2)-syr_true(2))^2);
    syr_weight(i) = (1 / sqrt(R) / sqrt(2*pi)) * exp(-d^2 / 2 / R); 
end
wsum = sum(syr_weight);
for i = 1:N
    syr_weight(i) = syr_weight(i)/wsum;
end
    utmp = rand;
for i = 1:N
    u = ((i-1) + utmp)/N;
    xtemsum = 0;
    for j = 1:N
        xtemsum = xtemsum + syr_weight(j);
        if xtemsum >=u
            syr_xpart(i,1) = syr_xarr(j,1);
            syr_xpart(i,2) = syr_xarr(j,2);
            syrw(i) = syr_weight(i);
            k = k+1;
            break;
        end
    end
end
%%
%更新结果
syr_xx = mean(syr_xpart);