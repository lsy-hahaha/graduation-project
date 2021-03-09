function [Hhat_EVD] = functionChannelEstimateEVD(H,nbrOfRealizations,M,K,L,p,f,tau_p,N)

%Generate pilot pattern  L是基站的数目 生成导频
if f == 1
    
    pilotPattern = ones(L,1);    %pilotPattern=（L×1）全1阵
    
elseif f == 2 %Only works in the running example with its 16 BSs
    
    pilotPattern = kron(ones(2,1),[1; 2; 1; 2; 2; 1; 2; 1]); %pilotPattern=（1 2 1 2 2 1 2 1 1 2 1 2 2 1 2 1）'
    
elseif f == 4 %Only works in the running example with its 16 BSs
    
    pilotPattern = kron(ones(2,1),[1; 2; 1; 2; 3; 4; 3; 4]); %pilotPattern=（1 2 1 2 3 4 3 4 1 2 1 2 3 4 3 4）'
    
elseif f == 16 %Only works in the running example with its 16 BSs
    
    pilotPattern = (1:L)'; %pilotPattern =(1 2 3 4 5 6……L)'
    
end



%Generate realizations of normalized noise 生成归一化噪声
Np = sqrt(0.5)*(randn(M,nbrOfRealizations,K,L,f) + 1i*randn(M,nbrOfRealizations,K,L,f));


%Prepare to store EVD channel estimates 用于存储EVD信道估计
Hhat_EVD= zeros(M,nbrOfRealizations,K,L);

% Go through all cells 遍历小区数目
for j = 1:L
    %Go through all f pilot groups 遍历所有导频组
    for g = 1:f
        
        %Extract the cells that belong to pilot group g 提取属于导频组g的小区
        %这块要是报错就取消注释
        %groupMembers = find(g==pilotPattern)';
        
        %Compute processed pilot signal for all UEs that use these pilots,
        %according to (5) 处理后的导频
        ypx = sqrt(p)*tau_p*sum(H(:,:,:,g==pilotPattern,j),4) + sqrt(tau_p)*Np(:,:,:,j,g);
    end
end 
%采样数N=50
%PLAN A对nbrOfRealizations进行采样 先取前50
Rypx=[];
ypx_sample=ypx(:,1:50,:,:,:);
for n = 1:N
Rypx1=(ypx_sample(:,n,:,:,:))*(ypx_sample(:,n,:,:,:))';
Rypx=Rypx1+Rypx;
end
R_y=(1/N)*Rypx; %把Ry尖求出来
%对Ry尖进行EVD分解
for i = 1:K
    for j = 1:L
        for k = 1:L
            [~,U_N(i,j,k)]=eig(R_y);
        end
    end
end

%定义y_n
%首先求矩阵的实部和虚部ypx_r和ypx_i
%把需要定义的函数直接写进来吧
for n = 1:N
ypx_new(:,n,:,:,:)=(ypx(:,n,:,:,:).')'; %ypx_new是ypx的共轭
ypx_r(:,n,:,:,:)=0.5*(ypx_new(:,n,:,:,:)+ypx(:,n,:,:,:));%实部矩阵
ypx_i(:,n,:,:,:)=0.5*(ypx(:,n,:,:,:)-ypx_new(:,n,:,:,:));%虚部矩阵
y_n(:,n,:,:,:)=[(ypx_r).' (ypx_i).']; %按列取数 相当于实部/虚部矩阵仅转置之后放在一起成为一个大矩阵
end
%定义变量An拔
%调用函数得到大尺度衰落系数矩阵
R = functionExampleSetup(L,K,M,ASDdeg,allLoS); 
%定义Sn拔
S_n=sqrt(p)*tau_p;
%定义An
for i = 1:K
    for j = 1:L
        for k=1:L
            beta(i,j,k)=trace(R(:,:,i,j,k));
        end
    end
end
An = sqrt(p)*U_N*R*beta*S_n;
%定义A_n拔 求An的实部虚部
A=[];
A_out=[];
for k = 1:L
An_r(:,:,k) = 0.5*((((An(:,:,k))').')+An);  
An_i(:,:,k) = 0.5*((((An(:,:,k))').')-An);  
A_n(:,:,k)  = [An_r(:,:,k) -An_i(:,:,k) ; An_i(:,:,k) An_r(:,:,k)];
%求ξ尖拔
%先把逆里面的求和写出来吧
A=A+A_n(:,:,k).'*A_n(:,:,k)
%再把伪逆外面的求和部分写出来
A_out =A_out+A_n(:,:,k).'*y_n(:,:,:,:,k);
xi_jian= (pinv(A))*A_out;
%表示单位矩阵
I_k=eye(K);
xijian=[I_k 1i*I_k]*xi_jian;
%求出模糊因子theta
theta_1=(diag(xijian)).';
theta=diag(theta,0)
