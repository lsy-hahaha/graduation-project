function [Hhat_EVD] = functionChannelEstimateEVD(H,nbrOfRealizations,M,K,L,p,f,tau_p,N)

%Generate pilot pattern  L�ǻ�վ����Ŀ ���ɵ�Ƶ
if f == 1
    
    pilotPattern = ones(L,1);    %pilotPattern=��L��1��ȫ1��
    
elseif f == 2 %Only works in the running example with its 16 BSs
    
    pilotPattern = kron(ones(2,1),[1; 2; 1; 2; 2; 1; 2; 1]); %pilotPattern=��1 2 1 2 2 1 2 1 1 2 1 2 2 1 2 1��'
    
elseif f == 4 %Only works in the running example with its 16 BSs
    
    pilotPattern = kron(ones(2,1),[1; 2; 1; 2; 3; 4; 3; 4]); %pilotPattern=��1 2 1 2 3 4 3 4 1 2 1 2 3 4 3 4��'
    
elseif f == 16 %Only works in the running example with its 16 BSs
    
    pilotPattern = (1:L)'; %pilotPattern =(1 2 3 4 5 6����L)'
    
end



%Generate realizations of normalized noise ���ɹ�һ������
Np = sqrt(0.5)*(randn(M,nbrOfRealizations,K,L,f) + 1i*randn(M,nbrOfRealizations,K,L,f));


%Prepare to store EVD channel estimates ���ڴ洢EVD�ŵ�����
Hhat_EVD= zeros(M,nbrOfRealizations,K,L);

% Go through all cells ����С����Ŀ
for j = 1:L
    %Go through all f pilot groups �������е�Ƶ��
    for g = 1:f
        
        %Extract the cells that belong to pilot group g ��ȡ���ڵ�Ƶ��g��С��
        %���Ҫ�Ǳ����ȡ��ע��
        %groupMembers = find(g==pilotPattern)';
        
        %Compute processed pilot signal for all UEs that use these pilots,
        %according to (5) �����ĵ�Ƶ
        ypx = sqrt(p)*tau_p*sum(H(:,:,:,g==pilotPattern,j),4) + sqrt(tau_p)*Np(:,:,:,j,g);
    end
end 
%������N=50
%PLAN A��nbrOfRealizations���в��� ��ȡǰ50
Rypx=[];
ypx_sample=ypx(:,1:50,:,:,:);
for n = 1:N
Rypx1=(ypx_sample(:,n,:,:,:))*(ypx_sample(:,n,:,:,:))';
Rypx=Rypx1+Rypx;
end
R_y=(1/N)*Rypx; %��Ry�������
%��Ry�����EVD�ֽ�
for i = 1:K
    for j = 1:L
        for k = 1:L
            [~,U_N(i,j,k)]=eig(R_y);
        end
    end
end

%����y_n
%����������ʵ�����鲿ypx_r��ypx_i
%����Ҫ����ĺ���ֱ��д������
for n = 1:N
ypx_new(:,n,:,:,:)=(ypx(:,n,:,:,:).')'; %ypx_new��ypx�Ĺ���
ypx_r(:,n,:,:,:)=0.5*(ypx_new(:,n,:,:,:)+ypx(:,n,:,:,:));%ʵ������
ypx_i(:,n,:,:,:)=0.5*(ypx(:,n,:,:,:)-ypx_new(:,n,:,:,:));%�鲿����
y_n(:,n,:,:,:)=[(ypx_r).' (ypx_i).']; %����ȡ�� �൱��ʵ��/�鲿�����ת��֮�����һ���Ϊһ�������
end
%�������An��
%���ú����õ���߶�˥��ϵ������
R = functionExampleSetup(L,K,M,ASDdeg,allLoS); 
%����Sn��
S_n=sqrt(p)*tau_p;
%����An
for i = 1:K
    for j = 1:L
        for k=1:L
            beta(i,j,k)=trace(R(:,:,i,j,k));
        end
    end
end
An = sqrt(p)*U_N*R*beta*S_n;
%����A_n�� ��An��ʵ���鲿
A=[];
A_out=[];
for k = 1:L
An_r(:,:,k) = 0.5*((((An(:,:,k))').')+An);  
An_i(:,:,k) = 0.5*((((An(:,:,k))').')-An);  
A_n(:,:,k)  = [An_r(:,:,k) -An_i(:,:,k) ; An_i(:,:,k) An_r(:,:,k)];
%��μ��
%�Ȱ�����������д������
A=A+A_n(:,:,k).'*A_n(:,:,k)
%�ٰ�α���������Ͳ���д����
A_out =A_out+A_n(:,:,k).'*y_n(:,:,:,:,k);
xi_jian= (pinv(A))*A_out;
%��ʾ��λ����
I_k=eye(K);
xijian=[I_k 1i*I_k]*xi_jian;
%���ģ������theta
theta_1=(diag(xijian)).';
theta=diag(theta,0)
