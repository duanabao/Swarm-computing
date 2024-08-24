function [X,Nt,Nh,count,xx,xpart,xpartw,xpartwr] = fun_resample(xarr,true,N,xpartw)
xpart = zeros(N,2);
X = zeros(N,2);
count = zeros(N,1);
Nh = 0;
Nl = 0; 
N1 = 0;
M = 0;
n = 1;
t = 0; weightMs = 0;
%%
%����Э�������
Covariance_matrix = cov(xarr);
[V, D] = eig(Covariance_matrix);
[largest_eigenvec_ind_c, r1] = find(D == max(max(D)));
%�����������
largest_c = V(:, largest_eigenvec_ind_c);
%�������ֵ
largest_v = max(max(D));
%��С����ֵ
if(largest_eigenvec_ind_c == 1)
    smallest_v = max(D(:,2));
else
    smallest_v = max(D(:,1));
end
%��Բ��б�Ƕ�
angle = atan2(largest_c(2), largest_c(1));
%���Ŷ�Ϊr1��r2,���ɶ�Ϊ2ʱ��Բ�Ĺ�ģ
r1=chi2inv(0.95,1);%0.85
r2=chi2inv(0.125,2);%0.1
% true = mean(xarr); 
%% �ز���
%Ȩ�ؼ��� %����1/N 
for i = 1:N
    d = ((xarr(i,1) - true(1))*cos(angle) + (xarr(i,2) - true(2))*sin(angle))^2/(largest_v) + ((xarr(i,2) - true(2))*cos(angle) - (xarr(i,1) - true(1))*sin(angle))^2/(smallest_v);
    if d > r1 
        Nl = Nl + 1;
        N1 = N1 + 1;
        count(N1) = i;
    elseif d >= r2 && d <= r1
        M = M + 1;
        xpart(i,1) = xarr(i,1);
        xpart(i,2) = xarr(i,2);
        weightMs = weightMs + xpartw(i);
        xpartw(i) = xpartw(i);
    else
        Nh = Nh + 1;
        X(Nh,1) = xarr(i,1);
        xpart(i,1) = xarr(i,1);
        X(Nh,2) = xarr(i,2);
        xpart(i,2) = xarr(i,2);
    end
end
weight = (1-weightMs)/(Nl+Nh);
% disp(Nl);
% disp(M);
% disp(Nh);
%%
%��������
Nt = Nl - floor(Nl/Nh)*Nh;
for i = 1:Nh   %�������С������
    if i <= Nt
        copytimes = floor(Nl/Nh) + 1;
        for j = 1:copytimes
            xpart(count(N1),1) = X(i,1);
            xpart(count(N1),2) = X(i,2);
            xpartw(count(N1)) = weight;
%             t = t+1;
%             Xtmp(t,1) = X(i,1);
%             Xtmp(t,2) = X(i,2);
            N1 = N1 - 1;
        end
    else
        copytimes = floor(Nl/Nh) + 1;    
        if copytimes ~=  1
            for j = 1:(copytimes-1)
                xpart(count(N1),1) = X(i,1);
                xpart(count(N1),2) = X(i,2);
                xpartw(count(N1)) = weight;
%                 t = t+1;
%                 Xtmp(t,1) = X(i,1);
%                 Xtmp(t,2) = X(i,2);
                N1 = N1 -1;
            end
        end
    end
end
xpartwr = xpartw;
%%
%����
xx = mean(xpart);