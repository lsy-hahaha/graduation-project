%Empty workspace and close figures
close all;
clear;
clc

%BS数目
L = 16;

%每个BS提供链路的UE数目
K = 5;

%BS天线数目
Mrange = 200;

%提取BS天线数目的最大值
Mmax = max(Mrange);

%定义导频复用因子
fRange = 1;

%选择随机终端位置设置的数量
nbrOfSetups = 7 ;

%选择每个设置能够实现的通道数
nbrOfRealizations = 100;

%传输带宽
B = 20e6;

%Total uplink transmit power per UE (W)
p = 0.1;
%BS处的噪声数值/dB
noiseFigure = 7;

%噪声方差/dBm
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;
%sigma2=db2pow(noiseVariancedBm-30); %W

%选择相关块的长度，由载波频率和外部因素决定
tau_c = 200;
%Angular standard deviation per path in the local scattering model (in degrees)
ASDdeg = 5;
%Power control parameter delta
deltadB=10;
%导频复用参数
f=1;
%导频长度
tau_p=K;
%Prelog factor assuming only UL transmission
prelogFactor=(tau_c -tau_p)/tau_c;


%存储仿真结果
userSE_MMSE= zeros(K*L,nbrOfSetups,1);
userSE_EWMMSE= zeros(K*L,nbrOfSetups,1);
userSE_LS= zeros(K*L,nbrOfSetups,1);

%Preallocate matrices for storing the covariance matrices and mean vectors
RNormalized=zeros(Mmax,Mmax,K,L,L);
HMeanNormalized=zeros(Mmax,K,L,L);


%% Go through all setups
for n = 1:nbrOfSetups
    
    %R and HMean is normalized i.e., norm(HMeanNormalized(:,k,l,j))^2=Mmax
    [RNormalized,HMeanNormalized,channelGaindB,ricianFactor,probLOS] = functionExampleSetup(L,K,Mmax,ASDdeg,0);
    
    %Controlled UL channel gain over noise
    channelGainUL= functionPowerControl( channelGaindB,noiseVariancedBm,deltadB,L);
    
    %If all UEs have same power
    % channelGainOverNoise=channelGaindB-noiseVariancedBm;
    
    for m=1:length(Mrange)
        
        %Generate the channel realizations and scale the normalized R and HMean
        [R_UL,HMean_UL,H_UL,H_UL_Rayleigh] = functionChannelGeneration( RNormalized(1:Mrange(m),1:Mrange(m),:,:,:)...
            ,HMeanNormalized(1:Mrange(m),:,:,:),channelGainUL,ricianFactor,probLOS,K,L,Mrange(m),nbrOfRealizations);
        
        %找到信道的理论矩阵
        Hhat_theoretical=H_UL;
        
        %The UL SE with MMSE for Rician fading
        %Generate the MMSE channel estimates
        Hhat_MMSE_UL = functionChannelEstimateMMSE(R_UL,HMean_UL,H_UL,nbrOfRealizations,Mrange(m),K,L,p,f,tau_p);
        %Compute UL SE using use and then forget bound in (19) with MMSE
        %SE_MonteCarlo_MMSE = functionMonteCarloSE_UL(Hhat_MMSE_UL,H_UL,prelogFactor,nbrOfRealizations,Mrange(m),K,L,p);
        
        
        
        %The UL SE with EW-MMSE for Rician fading
        %Generate the EWMMSE channel estimates
        Hhat_EWMMSE = functionChannelEstimateEWMMSE(R_UL,HMean_UL,H_UL,nbrOfRealizations,Mrange(m),K,L,p,f,tau_p);
        %Compute UL SE using use and then forget bound in (19) with EWMMSE
        %SE_MonteCarlo_EWMMSE= functionMonteCarloSE_UL(Hhat_EWMMSE,H_UL,prelogFactor,nbrOfRealizations,Mrange(m),K,L,p);
        
        
        
        %The UL SE with LS for Rician fading
        %Generate the LS channel estimates
        Hhat_LS = functionChannelEstimateLS(H_UL,nbrOfRealizations,Mrange(m),K,L,p,f,tau_p);
        %Compute UL SE using use and then forget bound in (19) with LS
        %SE_MonteCarlo_LS = functionMonteCarloSE_UL(Hhat_LS,H_UL,prelogFactor,nbrOfRealizations,Mrange(m),K,L,p);
        
        
        %Store SEs for all UEs
%        userSE_MMSE(:,n) = SE_MonteCarlo_MMSE(:);
%        userSE_EWMMSE(:,n) = SE_MonteCarlo_EWMMSE(:);
%        userSE_LS(:,n) = SE_MonteCarlo_LS(:);
     
        
        %Output simulation progress
        disp([num2str(Mrange(m)) ' antennas of ' num2str(Mmax)]);
        clear R_UL HMean_UL
    end
    
    %Output simulation progress
    disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
    clear RNormalized HMeanNormalized
    
end

%% Plot the figure
% CDFnumbers = linspace(0,1,K*L*nbrOfSetups);
% figure(1);
% hold on; box on;
% plot(sort(userSE_MMSE(:)),CDFnumbers,'r-');
% plot(sort(userSE_EWMMSE(:)),CDFnumbers,'k--');
% plot(sort(userSE_LS(:)),CDFnumbers,'b-.');
% xlabel('SE per UE [bit/s/Hz]');
% ylabel('Cumulative Distribution Function (CDF)');   %累积分布函数
% legend('MMSE estimator','EW-MMSE estimator','LS estimator','Location','SouthEast');

% H1=Hhat_MMSE_UL(:,:,1,1);
% H2=Hhat_EWMMSE(:,:,1,1);
NMSE1=[];
for y=1:Mrange
    for r=1:K
        for t=1:L
        NMSEt=norm(Hhat_LS(1:y,:,r,1,t)-H_UL(1:y,:,r,1,t));
        NMSE1(r,t,y)=NMSEt;
        end
    end
NMSE_w=mean2(NMSE1(:,:,y));
NMSE_f(y)=(mean(NMSE_w))/(100*y);
end
% H1=Hhat_LS(:,:,2,1);
% H_t=H_UL(:,:,2,1);

% NMSE1=mean(norm(H1-H_t));
% NMSE2=mean(norm(H2-H_t));


% for w=1:50
%     NMSE_w=mean(NMSE2);
%     NMSE=[NMSE NMSE_w];
% end
% NMSEx=NMSE(1,51:100);
% NMSE1=mean(norm(H3(1,:)-H_t(1,:)));
% NMSE2=mean(norm(H3(1:2,:)-H_t(1:2,:)));
% NMSE3=mean(norm(H3(1:3,:)-H_t(1:3,:)));
% NMSE4=mean(norm(H3(1:4,:)-H_t(1:4,:)));
% NMSE5=mean(norm(H3(1:5,:)-H_t(1:5,:)));
% NMSE6=mean(norm(H3(1:6,:)-H_t(1:6,:)));
% NMSE7=mean(norm(H3(1:7,:)-H_t(1:7,:)));
% NMSE8=mean(norm(H3(1:8,:)-H_t(1:8,:)));
% NMSE9=mean(norm(H3(1:9,:)-H_t(1:9,:)));
% NMSE10=mean(norm(H3-H_t));
% NMSE=[NMSE1 NMSE2 NMSE3 NMSE4 NMSE5,NMSE6 NMSE7 NMSE8 NMSE9 NMSE10];

  figure(1)
% % hold on; box on;
% % plot(1:size(H1,1),NMSE1,'r-');
% % plot(1:size(H2,1),NMSE2,'k--');
 plot(1:size(NMSE_f,2),NMSE_f,'r-');
 axis([100 200 0 0.07]);
% 
 xlabel('天线数目');
 ylabel('MSE'); 
% %legend('MMSE','EW-MMSE','LS','Location','SouthEast');