function [Hhat_LS] = functionChannelEstimateLS(H,nbrOfRealizations,M,K,L,p,f,tau_p)
%Generating LS channel estimates
%
%This Matlab function was developed to generate simulation results to:
%
%Ozgecan Ozdogan, Emil Bjornson, Erik G. Larsson, 鈥淢assive MIMO with
%Spatially Correlated Rician Fading Channels,鈥?IEEE Transactions on
%Communications, To appear.
%
%Download article: https://arxiv.org/abs/1805.07972
%
%This is version 1.0 (Last edited: 2019-02-01)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.


%Generate pilot pattern  L是基站的数目
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


%Prepare to store LS channel estimates LS
Hhat_LS= zeros(M,nbrOfRealizations,K,L);

% Go through all cells 遍历小区数目
for j = 1:L
    %Go through all f pilot groups 遍历所有导频组
    for g = 1:f
        
        %Extract the cells that belong to pilot group g 提取属于导频组g的小区
        groupMembers = find(g==pilotPattern)';
        
        %Compute processed pilot signal for all UEs that use these pilots,
        %according to (5) 处理后的导频
        ypx = sqrt(p)*tau_p*sum(H(:,:,:,g==pilotPattern,j),4) + sqrt(tau_p)*Np(:,:,:,j,g);
        %Go through all UEs 遍历所有用户
        for k=1:K
            %Go through the cells in pilot group g
            for l = groupMembers
                %Compute (15)
                Hhat_LS(:,:,k,l,j)= (1/(sqrt(p)*tau_p))*ypx(:,:,k);
            end
        end
    end
end
