function [Y_out,sme2] = G_CMAC( Sample, Output, Testsp, TestOut )
%[测试数据集的测试结果，测试的平均误差] = G_CMAC(训练数据集，训练数据输出集，测试数据集，测试数据输出集)
%数据集都是以列代表一个样本

xite=0.5;               %y和ym的偏差系数带来的权系数w的更新的比例系数η，参见42页式2.16
alfa=0.05;              %权系数的前一步和前两步的变化差带来的比例系数，这个书上没有
trainnum = 200;         %训练次数
pfErr = zeros(size(trainnum)); %%用于记录收敛情况

% N=759375;                    %权相量的维数，是不是就是A到AP映射的维数？m*nb^4
N = 14741;       %暂未用到，不必管
num = 16000;     %暂未用到，不必管
m = 8;                 %级数
nb = 7;                 %每级包含的块数 
M=(m * (nb - 1)+1);                  %量化时用到的系数
Wnum = m*nb^size(Sample,1);         %权值矩阵的大小
uc = 0.5;                   %CMAC中每块边界值在高斯函数中的值，用于计算高斯标准差
% gv = 9;                  %高斯标准差
% gv = m/(2*sqrt(-log(uc)));
xmin = -1;          %输入的最小最大值，用于量化计算
xmax = 1;           %
% w0=zeros(size(Output,1),N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%高斯函数期望m设置%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gu = zeros(m,nb,size(Sample,1));   %CMAC每块高斯函数期望值
Gv = ones(m,nb,size(Sample,1)).*m/(2*sqrt(-log(uc)));  %标准差。根据公式计算

%初始化高斯期望，取块的中间值为初始期望
for rm = 1:1:m
   for rn = 1:1:nb
       Gu(rm,rn,:) = (m/2-m+rm) + (rn-1)*m;
   end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w0 = zeros(size(Output,1),Wnum);  %权值矩阵
preci=1;  %记录训练次数
w=w0;
w_1=w0;
w_2=w0;
err = zeros(1,size(Output,1));  %误差，用于更新权值等参数
%err0 = 0;
ym = zeros(size(Output,1),size(Sample,2));   %训练样本经CMAC后的输出
Y_out = zeros(size(Output,1),size(Testsp,2));  %测试样本经CMAC后的输出
while (preci <= trainnum)  
  sme0 = 0;  %用于保存平方误差
 for k = 1:1:size(Sample,2)   %遍历每个训练样本
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%训练CMAC网络%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%注意这里的输入中有上一次计算出的权向量w%%%%%%%%%%%%%%%%%%
    ym(:,k)=G_CMACout(Sample(:,k),M,N,m,nb,w_1,Gu,Gv,xmin,xmax); %样本经CMAC后的输出
    
    for i = 1:1:size(Output,1)  %计算输出的每一位的误差。或许不用for循环就可以，未来得及优化
        err(i) = Output(i,k)-ym(i,k);                        
    end
    
    sme0 = sme0 + err*err'; %平方误差累加
    
    %G_CMACupdate为更新权值、高斯函数等参数
    [w,Gu,Gv]=G_CMACupdate(Sample(:,k),err,M,N,m,nb,w_1,w_2,Gu,Gv,xite,alfa,xmin,xmax);
    w_2=w_1;                   %权向量的前两次值，进行更新
    w_1=w;                     %权向量的前一次值，进行更新
 end
 
sme0 = sqrt(sme0/size(Sample,2));   %平均平方误差
 pfErr(preci) = sme0;  %记录每一步

 preci =preci +1;
end
epochs = linspace(1,trainnum,trainnum);%步数
figure(4);
hold on;
plot(epochs,pfErr);%画收敛图

sme = 0; %暂未用到
sme1 = 0; %暂未用到
sme2 = 0; %测试数据的平均平方误差
for k1 = 1:1:size(Testsp,2)
   Y_out(:,k1) = G_CMACout(Testsp(:,k1),M,N,m,nb,w,Gu,Gv,xmin,xmax); 
   e = 0;
   for i = 1:1:size(Output,1)
      e =  e + (TestOut(i,k1) - Y_out(i,k1))^2;
   end
   sme = sme + sqrt(e)/size(Output,1);%暂未用到
   sme1 = sme1 + sqrt(e);%暂未用到
   sme2 = sme2 + e;
end
sme = sme/size(Testsp,2);%暂未用到
sme1 = sme1/size(Testsp,2);%暂未用到
sme2 = sqrt(sme2/size(Testsp,2));%测试数据的平均平方误差
w15 = w;
gu = Gu;
gv = Gv;
%%保存CMAC的参数
save me_w w15;
save me_Gu gu;
save me_Gv gv;