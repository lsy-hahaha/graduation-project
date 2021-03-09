function [R,HMean,H,H_Rayleigh] = functionChannelGeneration( RNormalized,HMeanNormalized,channelGaindB,ricianFactor,probLOS,K,L,M,nbrOfRealizations)

%Prepare to store
channelGain_LOS=zeros(K,L,L);
channelGain_NLOS=zeros(K,L,L);
HMean=zeros(M,K,L,L);
R=zeros(M,M,K,L,L);

%Go through all UEs and apply the channel gains to the spatial
%correlation matrices and mean vectors
for k=1:K
    for l=1:L
        for j=1:L
            
            if probLOS(k,l,j)==1 %The LoS Path exists, Rician Factor ~= 0 ����LOS
                channelGain_LOS(k,l,j)= sqrt(ricianFactor(k,l,j)/(ricianFactor(k,l,j) +1 ))*db2pow(channelGaindB(k,l,j));
                channelGain_NLOS(k,l,j)=sqrt(1/(ricianFactor(k,l,j) +1 ))*db2pow(channelGaindB(k,l,j));
            else  %Pure NLoS case ��NLOS
                channelGain_LOS(k,l,j)= 0;
                channelGain_NLOS(k,l,j)=db2pow(channelGaindB(k,l,j));
            end
            %Scaling operation
            HMean(:,k,l,j)=sqrt(channelGain_LOS(k,l,j))*HMeanNormalized(:,k,l,j);
            R(:,:,k,l,j)=channelGain_NLOS(k,l,j)*RNormalized(:,:,k,l,j);
            
        end
    end
end


%Generate the channel realizations
%Generate uncorrelated random variables ���ɲ���ص��������
W = (randn(M,nbrOfRealizations,K,L,L)+1i*randn(M,nbrOfRealizations,K,L,L));

%Prepare to store channel realizations ׼���洢
H=zeros(M,nbrOfRealizations,K,L,L);
H_Rayleigh=zeros(M,nbrOfRealizations,K,L,L);

%Reshape the mean vectors to get same mean vector for all realizations
HMeanx=reshape(repmat(HMean,nbrOfRealizations,1),M,nbrOfRealizations,K,L,L);

%Go through all UEs and generate channel realizations
for j = 1:L
    
    for l = 1:L
        
        for k = 1:K
            Rsqrt = sqrtm(R(:,:,k,j,l));
            H(:,:,k,j,l) = sqrt(0.5)*Rsqrt*W(:,:,k,j,l) + HMeanx(:,:,k,j,l);
            H_Rayleigh(:,:,k,j,l) = sqrt(0.5)*Rsqrt*W(:,:,k,j,l);
        end
        
    end
    
end

end
