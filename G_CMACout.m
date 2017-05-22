function Out = G_CMACout( u,M,N,m,n,w,gu,gv,xmin,xmax )
%G_CMACout为直接使用CMAC
%   
d = 7;%步长，暂未用到
Out = zeros(size(w,1),1);
intsize = size(u,1);%输入的维数
gos = ones(1,m);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for j = 1:1:m  %遍历CMAC的每层
    %    for k = 1:1:n
        ad(j) = j;  %用于计算权值地址
          for i = 1:1:intsize  %遍历输入的每维
             t(i)=(u(i)-xmin)*(M)/(xmax-xmin);  %量化
    %          gu = (m/2-m+j)+(k-1)*m;
             k(i) = floor((t(i) + m - j)/m);  %计算输入属于该层的哪一块
             %以下两个if防止溢出块数之外
             if k(i)<0
                 k(i) = 0;
             end
             if k(i)>=n
                 k(i) = n-1;%%%?????
             end
             %计算在该块中的高斯函数值，每个输入累乘，参考导弹追踪论文里的bk(x)
             gos(j) = gos(j)*exp(-(t(i)-gu(j,k(i)+1,i))^2/gv(j,k(i)+1,i)^2);
             ad(j) = ad(j) + m*k(i)*n^(intsize-i); %计算权值地址
          end
    %    end
%         ad(j) = k(1)*n^3*m + k(2)*n^2*m + k(3)*n*m + k(4)*m + j;
        
        for o = 1:1:size(w,1) %计算输出
            Out(o) = Out(o) + w(o,ad(j));
        end
    end
    
end


