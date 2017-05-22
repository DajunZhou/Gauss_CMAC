function [Y_out,sme2] = G_CMAC( Sample, Output, Testsp, TestOut )
%[�������ݼ��Ĳ��Խ�������Ե�ƽ�����] = G_CMAC(ѵ�����ݼ���ѵ��������������������ݼ����������������)
%���ݼ��������д���һ������

xite=0.5;               %y��ym��ƫ��ϵ��������Ȩϵ��w�ĸ��µı���ϵ���ǣ��μ�42ҳʽ2.16
alfa=0.05;              %Ȩϵ����ǰһ����ǰ�����ı仯������ı���ϵ�����������û��
trainnum = 200;         %ѵ������
pfErr = zeros(size(trainnum)); %%���ڼ�¼�������

% N=759375;                    %Ȩ������ά�����ǲ��Ǿ���A��APӳ���ά����m*nb^4
N = 14741;       %��δ�õ������ع�
num = 16000;     %��δ�õ������ع�
m = 8;                 %����
nb = 7;                 %ÿ�������Ŀ��� 
M=(m * (nb - 1)+1);                  %����ʱ�õ���ϵ��
Wnum = m*nb^size(Sample,1);         %Ȩֵ����Ĵ�С
uc = 0.5;                   %CMAC��ÿ��߽�ֵ�ڸ�˹�����е�ֵ�����ڼ����˹��׼��
% gv = 9;                  %��˹��׼��
% gv = m/(2*sqrt(-log(uc)));
xmin = -1;          %�������С���ֵ��������������
xmax = 1;           %
% w0=zeros(size(Output,1),N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��˹��������m����%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gu = zeros(m,nb,size(Sample,1));   %CMACÿ���˹��������ֵ
Gv = ones(m,nb,size(Sample,1)).*m/(2*sqrt(-log(uc)));  %��׼����ݹ�ʽ����

%��ʼ����˹������ȡ����м�ֵΪ��ʼ����
for rm = 1:1:m
   for rn = 1:1:nb
       Gu(rm,rn,:) = (m/2-m+rm) + (rn-1)*m;
   end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w0 = zeros(size(Output,1),Wnum);  %Ȩֵ����
preci=1;  %��¼ѵ������
w=w0;
w_1=w0;
w_2=w0;
err = zeros(1,size(Output,1));  %�����ڸ���Ȩֵ�Ȳ���
%err0 = 0;
ym = zeros(size(Output,1),size(Sample,2));   %ѵ��������CMAC������
Y_out = zeros(size(Output,1),size(Testsp,2));  %����������CMAC������
while (preci <= trainnum)  
  sme0 = 0;  %���ڱ���ƽ�����
 for k = 1:1:size(Sample,2)   %����ÿ��ѵ������
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ѵ��CMAC����%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%ע�����������������һ�μ������Ȩ����w%%%%%%%%%%%%%%%%%%
    ym(:,k)=G_CMACout(Sample(:,k),M,N,m,nb,w_1,Gu,Gv,xmin,xmax); %������CMAC������
    
    for i = 1:1:size(Output,1)  %���������ÿһλ����������forѭ���Ϳ��ԣ�δ���ü��Ż�
        err(i) = Output(i,k)-ym(i,k);                        
    end
    
    sme0 = sme0 + err*err'; %ƽ������ۼ�
    
    %G_CMACupdateΪ����Ȩֵ����˹�����Ȳ���
    [w,Gu,Gv]=G_CMACupdate(Sample(:,k),err,M,N,m,nb,w_1,w_2,Gu,Gv,xite,alfa,xmin,xmax);
    w_2=w_1;                   %Ȩ������ǰ����ֵ�����и���
    w_1=w;                     %Ȩ������ǰһ��ֵ�����и���
 end
 
sme0 = sqrt(sme0/size(Sample,2));   %ƽ��ƽ�����
 pfErr(preci) = sme0;  %��¼ÿһ��

 preci =preci +1;
end
epochs = linspace(1,trainnum,trainnum);%����
figure(4);
hold on;
plot(epochs,pfErr);%������ͼ

sme = 0; %��δ�õ�
sme1 = 0; %��δ�õ�
sme2 = 0; %�������ݵ�ƽ��ƽ�����
for k1 = 1:1:size(Testsp,2)
   Y_out(:,k1) = G_CMACout(Testsp(:,k1),M,N,m,nb,w,Gu,Gv,xmin,xmax); 
   e = 0;
   for i = 1:1:size(Output,1)
      e =  e + (TestOut(i,k1) - Y_out(i,k1))^2;
   end
   sme = sme + sqrt(e)/size(Output,1);%��δ�õ�
   sme1 = sme1 + sqrt(e);%��δ�õ�
   sme2 = sme2 + e;
end
sme = sme/size(Testsp,2);%��δ�õ�
sme1 = sme1/size(Testsp,2);%��δ�õ�
sme2 = sqrt(sme2/size(Testsp,2));%�������ݵ�ƽ��ƽ�����
w15 = w;
gu = Gu;
gv = Gv;
%%����CMAC�Ĳ���
save me_w w15;
save me_Gu gu;
save me_Gv gv;