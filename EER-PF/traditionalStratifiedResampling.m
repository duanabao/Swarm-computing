function [sr_xx,sr_xpart,srw] = traditionalStratifiedResampling(sr_xarr,sr_true,N)
%%
%��������
sr_xpart = zeros(N,2);%����������Ϣ
xweight = zeros(N,1);%����Ȩ��
R =0.1;
k =0;
%%
%�㷨ʵ��
for i = 1:N
    d = sqrt((sr_xarr(i,1) - sr_true(1))^2+(sr_xarr(i,2)-sr_true(2))^2);
    xweight(i) = (1 / sqrt(R) / sqrt(2*pi)) * exp(-d^2 / 2 / R); 
end
wsum = sum(xweight);
for i = 1:N
    xweight(i) = xweight(i)/wsum;
end
for i = 1:N
    utmp = rand;
    u = ((i-1) + utmp)/N;
    xtemsum = 0;
    for j = 1:N
        xtemsum = xtemsum + xweight(j);
        if xtemsum >=u
            sr_xpart(i,1) = sr_xarr(j,1);
            sr_xpart(i,2) = sr_xarr(j,2);
            srw(i) = xweight(i);
            k = k+1;
            break;
        end
    end
end
%%
%���½��
sr_xx = mean(sr_xpart);