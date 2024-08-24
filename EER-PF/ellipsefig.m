 %画置信椭圆

    function  h = ellipsefig(xc,P,r) %r=5.991
    
    [V, D] = eig(P); %求P的全部特征值构成对角阵D，特征向量构成V
    
    %计算最大特征值和最大特征向量
    
    [largest_eigenvec_ind_c, r1] = find(D == max(max(D)));  %largest_eigenvec_ind_c存入最大特征值所在的位置
    
    largest_c = V(:, largest_eigenvec_ind_c);  %将对应位置的特征向量存到largest_c
    
    largest_v = max(max(D));
    
    %计算最小特征值和最小特征向量
    
    if(largest_eigenvec_ind_c == 1)
        
        smallest_v = max(D(:,2));
        
    else
        
        smallest_v = max(D(:,1));
        
    end
    
    %最大特征向量和x轴之间的夹角
    
    angle = atan2(largest_c(2), largest_c(1));
    
    if(angle < 0)
        
        angle = angle + 2*pi;
        
    end
    
    if(length(xc)>2)  %判断矩阵xc的行列数最大值是否大于2
        
        avg1 = mean(xc);  %返回xc每列的均值
        
    else
        
        avg1=xc;
        
    end
    
    q=sqrt(r);
    
    t = linspace(0, 2*pi);
    
    X0=avg1(1);
    
    Y0=avg1(2);
    
    aa = q*sqrt(largest_v);
    
    bb = q*sqrt(smallest_v);
    
    xr=aa*cos(t);
    
    yr=bb*sin(t);
    
    %旋转矩阵
    R = [ cos(angle) -sin(angle); sin(angle) cos(angle) ];
    xy = R*[xr;yr];  % 旋转椭圆
    h = plot(xy(1,:)+X0,xy(2,:)+Y0, '-.','linewidth',2);