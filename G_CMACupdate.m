function [w,gu,gv] = G_CMACupdate( u,err,M,N,m,n,w_1,w_2,gu,gv,xite,alfa,xmin,xmax )
%G_CMACupdate为更新参数

d = 7;%步长，暂未用到
w=w_1;
intsize = size(u,1);%输入的维数

gos = ones(1,m); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:1:m  %遍历CMAC的每层
%    for k = 1:1:n
       ad(j) = j; %用于计算权值地址
      for i = 1:1:intsize    %遍历输入的每维
         t(i)=(u(i)-xmin)*(M)/(xmax-xmin);  %量化
         k(i) = floor((t(i) + m - j)/m);   %计算输入属于该层的哪一块
         %以下两个if防止溢出块数之外
         if k(i)<0
             k(i) = 0;
         end
         if k(i)>=n
             k(i) = n-1;
         end
         %计算在该块中的高斯函数值，每个输入累乘，参考导弹追踪论文里的bk(x)
         gos(j) = gos(j)*exp(-(t(i)-gu(j,k(i)+1,i))^2/gv(j,k(i)+1,i)^2);
         ad(j) = ad(j) + m*k(i)*n^(intsize-i); %计算权值地址
      end
%    end
    
    for o = 1:1:size(w,1)     %权值更新           
        d_w(o)=xite*err(o)/m;%*gos(j);
        w(o,ad(j)) = w(o,ad(j))+ d_w(o);   
    end
    
    %for，高斯函数更新
%     for i0 = 1:1:intsize
%         for o1 = 1:1:size(w,1)
%             nw(o1) = w(o1,ad(j));
%         end
%         new_gu = gos(j)*(xite/m)*err*nw'*(2*(t(i0)-gu(j,k(i0)+1,i0))/gv(j,k(i0)+1,i0)^2);%除以m
%         new_gv = gos(j)*(xite/m)*err*nw'*(2*(t(i0)-gu(j,k(i0)+1,i0))^2/gv(j,k(i0)+1,i0)^3);
%         gu(j,k(i0)+1,i0) = gu(j,k(i0)+1,i0) + new_gu;
%         gv(j,k(i0)+1,i0) = gv(j,k(i0)+1,i0) +new_gv;
%     end
end
end
