function Out = G_CMACout( u,M,N,m,n,w,gu,gv,xmin,xmax )
%G_CMACoutΪֱ��ʹ��CMAC
%   
d = 7;%��������δ�õ�
Out = zeros(size(w,1),1);
intsize = size(u,1);%�����ά��
gos = ones(1,m);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for j = 1:1:m  %����CMAC��ÿ��
    %    for k = 1:1:n
        ad(j) = j;  %���ڼ���Ȩֵ��ַ
          for i = 1:1:intsize  %���������ÿά
             t(i)=(u(i)-xmin)*(M)/(xmax-xmin);  %����
    %          gu = (m/2-m+j)+(k-1)*m;
             k(i) = floor((t(i) + m - j)/m);  %�����������ڸò����һ��
             %��������if��ֹ�������֮��
             if k(i)<0
                 k(i) = 0;
             end
             if k(i)>=n
                 k(i) = n-1;%%%?????
             end
             %�����ڸÿ��еĸ�˹����ֵ��ÿ�������۳ˣ��ο�����׷���������bk(x)
             gos(j) = gos(j)*exp(-(t(i)-gu(j,k(i)+1,i))^2/gv(j,k(i)+1,i)^2);
             ad(j) = ad(j) + m*k(i)*n^(intsize-i); %����Ȩֵ��ַ
          end
    %    end
%         ad(j) = k(1)*n^3*m + k(2)*n^2*m + k(3)*n*m + k(4)*m + j;
        
        for o = 1:1:size(w,1) %�������
            Out(o) = Out(o) + w(o,ad(j));
        end
    end
    
end


