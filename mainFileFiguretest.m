%Empty workspace and close figures
close all;
clear;
clc

%BS��Ŀ
L = 16;

%ÿ��BS�ṩ��·��UE��Ŀ
K = 5;

%BS������Ŀ
Mrange = 200;

%��ȡBS������Ŀ�����ֵ
Mmax = max(Mrange);

%���嵼Ƶ��������
fRange = 1;

%ѡ������ն�λ�����õ�����
nbrOfSetups = 7 ;

%ѡ��ÿ�������ܹ�ʵ�ֵ�ͨ����
nbrOfRealizations = 100;

%�������
B = 20e6;

%Total uplink transmit power per UE (W)
% p = 0.1;
%BS����������ֵ/dBm
noiseFigure = 7;

%��������/dBm
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;
sigma2=db2pow(noiseVariancedBm-30); %W

%ѡ����ؿ�ĳ��ȣ����ز�Ƶ�ʺ��ⲿ���ؾ���
tau_c = 200;
%Angular standard deviation per path in the local scattering model (in degrees)
ASDdeg = 5;
%Power control parameter delta
deltadB=10;
%��Ƶ���ò���
f=2;
%��Ƶ����
tau_p=K;
%Prelog factor assuming only UL transmission
prelogFactor=(tau_c -tau_p)/tau_c;

%%�����Ƿ���������
%%�����Ǵ洢������

%�洢������
% userSE_MMSE= zeros(K*L,nbrOfSetups,1);
% userSE_EWMMSE= zeros(K*L,nbrOfSetups,1);
% userSE_LS= zeros(K*L,nbrOfSetups,1);

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
        
        %�ҵ��ŵ������۾���
        Hhat_theoretical=H_UL;
        
        for p = 1:5
        %The UL SE with MMSE for Rician fading
        %Generate the MMSE channel estimates
%         Hhat_MMSE_UL = functionChannelEstimateMMSE(R_UL,HMean_UL,H_UL,nbrOfRealizations,Mrange(m),K,L,p,f,tau_p);
        %Compute UL SE using use and then forget bound in (19) with MMSE
        %SE_MonteCarlo_MMSE = functionMonteCarloSE_UL(Hhat_MMSE_UL,H_UL,prelogFactor,nbrOfRealizations,Mrange(m),K,L,p);
        
        
        
        %The UL SE with EW-MMSE for Rician fading
        %Generate the EWMMSE channel estimates
%         Hhat_EWMMSE = functionChannelEstimateEWMMSE(R_UL,HMean_UL,H_UL,nbrOfRealizations,Mrange(m),K,L,p,f,tau_p);
        %Compute UL SE using use and then forget bound in (19) with EWMMSE
        %SE_MonteCarlo_EWMMSE= functionMonteCarloSE_UL(Hhat_EWMMSE,H_UL,prelogFactor,nbrOfRealizations,Mrange(m),K,L,p);
        
        
        
        %The UL SE with LS for Rician fading
        %Generate the LS channel estimates
        Hhat_LS = functionChannelEstimateLS(H_UL,nbrOfRealizations,Mrange(m),K,L,p,f,tau_p);
        
        NMSE1=[];
        for r=1:K
            for t=1:L
                NMSEt=norm(Hhat_LS(:,:,r,1,t)- Hhat_theoretical(:,:,r,1,t));
                 NMSE1(r,t)=NMSEt;
            end
        end
        NMSE(p)=mean2(NMSE1);
        end
                 
                
        

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
%     clear RNormalized HMeanNormalized
    
end




% NMSE1=[];
% for r= 1:K
%     for t=1:L
%         NMSEt=norm(Hhat_LS1(:,:,r,1,t)- Hhat_theoretical(:,:,r,1,t));
%         NMSE1(r,t)=NMSEt;
%     end
% end
%  NMSE=mean2(NMSE1);


%% Plot the figure
% CDFnumbers = linspace(0,1,K*L*nbrOfSetups);
% figure(1);
% hold on; box on;
% plot(sort(userSE_MMSE(:)),CDFnumbers,'r-');
% plot(sort(userSE_EWMMSE(:)),CDFnumbers,'k--');
% plot(sort(userSE_LS(:)),CDFnumbers,'b-.');
% xlabel('SE per UE [bit/s/Hz]');
% ylabel('Cumulative Distribution Function (CDF)');   %�ۻ��ֲ�����
% legend('MMSE estimator','EW-MMSE estimator','LS estimator','Location','SouthEast');





   figure(1)
% % % hold on; box on;
% % % plot(1:size(H1,1),NMSE1,'r-');
% % % plot(1:size(H2,1),NMSE2,'k--');
  plot([0.1 0.2 0.3 0.4 0.5] ,NMSE);
% 
  xlabel('p');
  ylabel('MSE'); 
% % %legend('MMSE','EW-MMSE','LS','Location','SouthEast');
% 

