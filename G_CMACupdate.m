function [w,gu,gv] = G_CMACupdate( u,err,M,N,m,n,w_1,w_2,gu,gv,xite,alfa,xmin,xmax )
%G_CMACupdateΪ���²���

d = 7;%��������δ�õ�
w=w_1;
intsize = size(u,1);%�����ά��

gos = ones(1,m); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:1:m  %����CMAC��ÿ��
%    for k = 1:1:n
       ad(j) = j; %���ڼ���Ȩֵ��ַ
      for i = 1:1:intsize    %���������ÿά
         t(i)=(u(i)-xmin)*(M)/(xmax-xmin);  %����
         k(i) = floor((t(i) + m - j)/m);   %�����������ڸò����һ��
         %��������if��ֹ�������֮��
         if k(i)<0
             k(i) = 0;
         end
         if k(i)>=n
             k(i) = n-1;
         end
         %�����ڸÿ��еĸ�˹����ֵ��ÿ�������۳ˣ��ο�����׷���������bk(x)
         gos(j) = gos(j)*exp(-(t(i)-gu(j,k(i)+1,i))^2/gv(j,k(i)+1,i)^2);
         ad(j) = ad(j) + m*k(i)*n^(intsize-i); %����Ȩֵ��ַ
      end
%    end
    
    for o = 1:1:size(w,1)     %Ȩֵ����           
        d_w(o)=xite*err(o)/m;%*gos(j);
        w(o,ad(j)) = w(o,ad(j))+ d_w(o);   
    end
    
    %for����˹��������
%     for i0 = 1:1:intsize
%         for o1 = 1:1:size(w,1)
%             nw(o1) = w(o1,ad(j));
%         end
%         new_gu = gos(j)*(xite/m)*err*nw'*(2*(t(i0)-gu(j,k(i0)+1,i0))/gv(j,k(i0)+1,i0)^2);%����m
%         new_gv = gos(j)*(xite/m)*err*nw'*(2*(t(i0)-gu(j,k(i0)+1,i0))^2/gv(j,k(i0)+1,i0)^3);
%         gu(j,k(i0)+1,i0) = gu(j,k(i0)+1,i0) + new_gu;
%         gv(j,k(i0)+1,i0) = gv(j,k(i0)+1,i0) +new_gv;
%     end
end
end
