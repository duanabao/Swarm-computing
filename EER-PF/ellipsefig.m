 %��������Բ

    function  h = ellipsefig(xc,P,r) %r=5.991
    
    [V, D] = eig(P); %��P��ȫ������ֵ���ɶԽ���D��������������V
    
    %�����������ֵ�������������
    
    [largest_eigenvec_ind_c, r1] = find(D == max(max(D)));  %largest_eigenvec_ind_c�����������ֵ���ڵ�λ��
    
    largest_c = V(:, largest_eigenvec_ind_c);  %����Ӧλ�õ����������浽largest_c
    
    largest_v = max(max(D));
    
    %������С����ֵ����С��������
    
    if(largest_eigenvec_ind_c == 1)
        
        smallest_v = max(D(:,2));
        
    else
        
        smallest_v = max(D(:,1));
        
    end
    
    %�������������x��֮��ļн�
    
    angle = atan2(largest_c(2), largest_c(1));
    
    if(angle < 0)
        
        angle = angle + 2*pi;
        
    end
    
    if(length(xc)>2)  %�жϾ���xc�����������ֵ�Ƿ����2
        
        avg1 = mean(xc);  %����xcÿ�еľ�ֵ
        
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
    
    %��ת����
    R = [ cos(angle) -sin(angle); sin(angle) cos(angle) ];
    xy = R*[xr;yr];  % ��ת��Բ
    h = plot(xy(1,:)+X0,xy(2,:)+Y0, '-.','linewidth',2);