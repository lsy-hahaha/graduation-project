function SE_Theoretical = functionTheoreticalSE_DL_MMSE( R,HMean,M,K,L,p,tau_p,prelogFactor,rho)
%Computes the DL SE with MMSE estimator in Theorem 4
%Note that the covariance matrices and the mean vectors are same for all
%channel realizations.
%
%This Matlab function was developed to generate simulation results to:
%
%Ozgecan Ozdogan, Emil Bjornson, Erik G. Larsson, “Massive MIMO with
%Spatially Correlated Rician Fading Channels,” IEEE Transactions on
%Communications, To appear.
%
%Download article: https://arxiv.org/abs/1805.07972
%
%This is version 1.0 (Last edited: 2019-02-01)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.


%Store identity matrix of size M x M
eyeM = eye(M);

%Prepare to store  (42) and (43)
CCterm1=zeros(K,L);
CCterm2=zeros(K,L);
CCterm2_p1=zeros(K,L,L,K);

%Prepare to store DL SE
SE_Theoretical=zeros(K,L);


% Go through all BSs
for j = 1:L
    
    %Go through all UEs
    for k = 1:K
        
        %Compute the matrix that is inverted in the MMSE estimator
        PsiInv = (p*tau_p*sum(R(:,:,k,:,j),4) + eyeM);
        %Compute (42) in Theorem 4
        CCterm1(k,j)=abs( p*tau_p*trace(R(:,:,k,j,j)/PsiInv*R(:,:,k,j,j)) + norm(HMean(:,k,j,j))^2  ) ;
        
        
        
        for lx=1:L
            for i=1:K
                
                %Compute the matrix that is inverted in the MMSE estimator for
                %all interfering UEs
                PsiInv_li = (p*tau_p*sum(R(:,:,i,:,lx),4) + eyeM);
                
                %Calculate (43) in Theorem 4
                %Compute the terms for non-pilot contaminated and pilot contaminated UEs
                CCterm2_p1(k,j,lx,i)=(p*tau_p*(trace(R(:,:,k,j,lx)*R(:,:,i,lx,lx)/PsiInv_li*R(:,:,i,lx,lx))) +...
                    (p*tau_p*HMean(:,k,j,lx)'*R(:,:,i,lx,lx)/PsiInv_li*R(:,:,i,lx,lx)*HMean(:,k,j,lx))+...
                    ( HMean(:,i,lx,lx)'*R(:,:,k,j,lx)*HMean(:,i,lx,lx))+...
                    abs(HMean(:,k,j,lx)'*HMean(:,i,lx,lx))^2)/...
                    abs( p*tau_p*trace(R(:,:,i,lx,lx)/PsiInv_li*R(:,:,i,lx,lx)) + norm(HMean(:,i,lx,lx))^2  );
                if k==i  %Compute the terms for pilot contaminated UEs
                    CCterm2_p1(k,j,lx,i)= CCterm2_p1(k,j,lx,i)+ ((p*p*tau_p*tau_p*abs(trace(R(:,:,k,j,lx)/PsiInv_li*R(:,:,i,lx,lx)))^2+...
                        2*p*tau_p*real(trace(R(:,:,k,j,lx)/PsiInv_li*R(:,:,i,lx,lx))*HMean(:,i,lx,lx)'*HMean(:,k,j,lx)))/...
                        abs( p*tau_p*trace(R(:,:,i,lx,lx)/PsiInv_li*R(:,:,i,lx,lx)) + norm(HMean(:,i,lx,lx))^2  ));
                end
                
            end
        end
        
        %Sum all the interferences from all UEs
        CCterm2(k,j)=(sum(sum(CCterm2_p1(k,j,:,:),3),4));
    end
end



%Compute the SE with (44) for each UE
for k=1:K
    for j=1:L
        SE_Theoretical(k,j)=prelogFactor*real(log2( 1 + (rho*CCterm1(k,j))/(rho*CCterm2(k,j)- rho*CCterm1(k,j) + 1 ))) ;
    end
end
