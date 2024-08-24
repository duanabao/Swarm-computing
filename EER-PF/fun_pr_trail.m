 function [xmen,xpr,sr_xpr,syr_xpr,xarr,xpartwr,srw,syrw] = fun_pr_trail(x,num,step,angel,step_mean,step_sigma,angel_mean,angel_sigma,N) 

global xpdisp xfdisp xcdisp Xtmp XNt XNh
xpart = zeros(N,2); %存储粒子
sr_xpart = zeros(N,2);
syr_xpart = zeros(N,2);

pr_step = zeros(num-1,1); %步长
sr_step = zeros(num-1,1);
syr_step = zeros(num-1,1);

pr_angel = zeros(num-1,1); %角度
sr_angel =zeros(num-1,1);
syr_angel =zeros(num-1,1);

xpr=zeros(num,2);  %预测轨迹点 
sr_xpr = zeros(num,2);
syr_xpr = zeros(num,2);

%第一点的位置生成
xpr(1,1) = x(1,1);
xpr(1,2) = x(1,2);
sr_xpr(1,1) = x(1,1);
sr_xpr(1,2) = x(1,2);
syr_xpr(1,1) = x(1,1);
syr_xpr(1,2) = x(1,2);

xarr = zeros(N,2); %存储粒子
sr_xarr = zeros(N,2);
syr_xarr = zeros(N,2);

%随机初始化N个粒子
for i = 1:N
    xpart(i,1) = x(1,1)+normrnd(0.1,0.1);
    xpart(i,2) = x(1,2)+normrnd(0.1,0.1);
end
for i = 1:N
    sr_xpart(i,1) = x(1,1)+normrnd(0.1,0.1);
    sr_xpart(i,2) = x(1,2)+normrnd(0.1,0.1);
end
for i = 1:N
    syr_xpart(i,1) = x(1,1)+normrnd(0.1,0.1);
    syr_xpart(i,2) = x(1,2)+normrnd(0.1,0.1);
    xpartw(1,i) = normrnd(10,0.5);
end
 qsum = sum(xpartw);
for i = 1 : N
    xpartw(i) = xpartw(i) / qsum;%权值归一化
end  
%粒子均值
xmeanpart = xpart;
%不使用重采样算法的第一个位置估计
xmen(1,1) = mean(xpart(:,1));
xmen(1,2) = mean(xpart(:,2));
%%
%使用本文提出的重采样算法对余下位置的预测,根据上一时刻的粒子状态决定
for i = 2:num
    tic;
    for j = 1:N
        pr_step(i-1) = step(i-1) + normrnd(step_mean,step_sigma);
        pr_angel(i-1) = rem(angel(i-1) + normrnd(angel_mean,angel_sigma),360);
        xarr(j,1) = xpart(j,1) + pr_step(i-1)*cos(pr_angel(i-1)/180*pi) + normrnd(0.1,0.1);
        xarr(j,2) = xpart(j,2) + pr_step(i-1)*sin(pr_angel(i-1)/180*pi) + normrnd(0.1,0.1);
    end
%     xarrtmp = xpart; %临时存储粒子数组
    true(1,1) = xpr(i-1,1)+pr_step(i-1)*cos(pr_angel(i-1)/180*pi) ; %选择置信椭圆中心点
    true(1,2) = xpr(i-1,2)+pr_step(i-1)*sin(pr_angel(i-1)/180*pi) ;
    if i == 6
        xfdisp = xarr;
    end
%     [X,Nt,Nh,count,xx,xpart] = fun_resample(xarr,true,N); %本文提出的PF重采样
    [X,Nt,Nh,count,xx,xpart,xpartw,xpartwr] = fun_resample(xarr,true,N,xpartw);
%     xpartw = xpartwr;
    xpr(i,1) = xx(1); %返回预测值
    xpr(i,2) = xx(2);
    if i == 6
        xpdisp = xpart;
        xcdisp = count;
        Xtmp = X;
        XNt = Nt;
        XNh = Nh;
    end
    toc;
    tEER(i) = toc;
end
%%
%使用SR重采样
for i = 2:num
    tic;
    for j = 1:N
        sr_step(i-1) = step(i-1) + normrnd(step_mean,step_sigma);
        sr_angel(i-1) = rem(angel(i-1) + normrnd(angel_mean,angel_sigma),360);
        sr_xarr(j,1) = sr_xpart(j,1) + sr_step(i-1)*cos(sr_angel(i-1)/180*pi) + normrnd(0.1,0.1);
        sr_xarr(j,2) = sr_xpart(j,2) + sr_step(i-1)*sin(sr_angel(i-1)/180*pi) + normrnd(0.1,0.1);
    end
    sr_true(1,1) = sr_xpr(i-1,1)+sr_step(i-1)*cos(sr_angel(i-1)/180*pi);
    sr_true(1,2) = sr_xpr(i-1,2)+sr_step(i-1)*sin(sr_angel(i-1)/180*pi);
    [sr_xx,sr_xpart,srw] = traditionalStratifiedResampling(sr_xarr,sr_true,N); %传统分层采样
    sr_xpr(i,1) = sr_xx(1); %返回预测值
    sr_xpr(i,2) = sr_xx(2);
    toc;
    tstr(i) = toc;
end
%%
%使用系统采样SR
for i = 2:num
    tic;
    for j = 1:N
        syr_step(i-1) = step(i-1) + normrnd(step_mean,step_sigma);
        syr_angel(i-1) = rem(angel(i-1) + normrnd(angel_mean,angel_sigma),360);
        syr_xarr(j,1) = syr_xpart(j,1) + syr_step(i-1)*cos(syr_angel(i-1)/180*pi) + normrnd(0.1,0.1);
        syr_xarr(j,2) = syr_xpart(j,2) + syr_step(i-1)*sin(syr_angel(i-1)/180*pi) + normrnd(0.1,0.1);
    end
    syr_true(1,1) = syr_xpr(i-1,1)+syr_step(i-1)*cos(syr_angel(i-1)/180*pi);
    syr_true(1,2) = syr_xpr(i-1,2)+syr_step(i-1)*sin(syr_angel(i-1)/180*pi);
    [syr_xx,syr_xpart,syrw] = SystematicResampling(syr_xarr,syr_true,N);
    syr_xpr(i,1) = syr_xx(1); %返回预测值
    syr_xpr(i,2) = syr_xx(2);
    toc;
    tsr(i) = toc;
end
%%
%不使用重采样算法的预测
for i = 2:num
    tic;
    for j = 1:N
        meanstep(i-1) = step(i-1) + normrnd(step_mean,step_sigma);
        meanangel(i-1) = rem(angel(i-1) + normrnd(angel_mean,angel_sigma),360);
        xmeanpart(j,1) = xmeanpart(j,1) + meanstep(i-1)*cos(meanangel(i-1)/180*pi) + normrnd(0.1,0.1);
        xmeanpart(j,2) = xmeanpart(j,2) + meanstep(i-1)*sin(meanangel(i-1)/180*pi) + normrnd(0.1,0.1);
    end
    xmen(i,1) = mean(xmeanpart(:,1));
    xmen(i,2) = mean(xmeanpart(:,2));
    toc;
    tur(i) = toc;
end
% save('count','count1','count2','count3');
% save('ResamplingTime','t12');
% save('time','tEER','tstr','tsr','tur');