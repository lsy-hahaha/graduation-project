clear;
close all;

%��С�����û�


Nd = 100;%���ų���

tao = 3;%��Ƶ����
r = 500;%С�����뾶
rc = r * 0.8;
rh = 100;%С����С�뾶
ra = rc / rh - 1;
gamma = 3.8;%·��˥��ָ��
mu = 0;%��Ӱ˥���ֵ
sigma = 10^0.8;%��Ӱ˥�䷽��
i_ant = 128;%���������Ŀ
Num = 100;%����
Nt = tao + Nd;%�ܳ���

SNR_n = 7;
SNR = [0;5;10;15;20;25;30];


H = zeros(i_ant,1);
G = zeros(i_ant,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mse_store = zeros(SNR_n, 3);

for ii = 1 : SNR_n
%     amp = 10 ^ (SNR(ii,1)*0.05) / sqrt(i_ant);
    bound = 0;
    bound1 = 0;
    csi_mse = zeros(1,3);
   for jj = 1 : Num
        shadow_amp = lognrnd(mu,sigma);%������Ӱ˥��

       %�����ŵ�
        H = 1 / sqrt(2) * (randn(i_ant,1)+1i*randn(i_ant,1));
        D =  shadow_amp  * ((rc*0.01).^(-0.5*gamma)); %��߶�˥�� ��ֵ
        G = H * D;

        beta = D; 
        amp_symbol = 10 ^ (SNR(ii,1)*0.05) / (sqrt(i_ant) * beta);
        amp_pilot = (sqrt(tao) * 10 ^ (SNR(ii,1)*0.05)) / (sqrt(i_ant) * beta);

        %���ɵ��Ʒ��ţ�ÿN�ж�Ӧһ��С���еķ���
        symbol = sign(randn(1,Nd));%BPSK����
            for l = 1 : Nd
                if 0 == symbol(1,l)
                    symbol(1,l) = 1;
                else
                end
            end

        receive_symbol = zeros(i_ant,Nd); %�洢�����ź�
        noise = (randn(i_ant,Nd)+1i*randn(i_ant,Nd)) / sqrt(2);%��������
        receive_symbol = receive_symbol + noise; %����
        receive_symbol(:,Nd) = receive_symbol(:,Nd) + amp_symbol * G * symbol(:,Nd);


        %% ���ɵ�Ƶ�ź� 
        pilots = zeros(tao,1);%Chu 
        if 0 == mod(tao,2)
            for k = 1 : 1
                for l = 0 : tao-1
                    ll = mod(l+k-1, tao);
                    pilots(l+1,k) = exp(1i*pi*ll*ll/tao);
                end
            end
        else
            for k = 1 : 1
                for l = 0 : tao-1
                    ll = mod(l+k-1, tao);
                    pilots(l+1,k) = exp(1i*pi*ll*(ll+1)/tao);
                end
            end
        end
        R_p = pilots.' * conj(pilots);
       %% ���յ�Ƶ�ź�
        receive_pilots = zeros(i_ant,tao);
        noise =  (randn(i_ant,tao)+1i*randn(i_ant,tao)) / sqrt(2);
        receive_pilots = receive_pilots + noise;
        receive_pilots = receive_pilots + amp_pilot * G  * pilots.';

        %% �ŵ�����
        H_svd = zeros(i_ant,1);
        H_ls = zeros(i_ant,1);
        theta = zeros(1,1);
        H_evd = zeros(i_ant,1);
        U = zeros(i_ant,1);
        
        r_x = zeros(i_ant,i_ant);
        for j = 1:tao
        r_x = r_x + receive_pilots(:,j) * receive_pilots(:,j)';
        end

        for k = 1 : Nd
            r_x = r_x + receive_symbol(:,k) * receive_symbol(:,k)';
        end
            r_x = r_x / Nt;   %������źŵ�Э������� ��Ƶ������
            %LS
            Hpi = receive_pilots * conj(pilots) / R_p / amp_symbol;
            H_ls = Hpi / D;
            
            %EVD 
            [U_all,~] = eig(r_x);
            [U_index,index] = sort(U_all,2);
            U(:,1) = U_index(:,128);
            A_sum = zeros(2,2);
            Ay = zeros(2,1);
            for n = 1 : tao
                yn = [real(receive_pilots(:,n));imag(receive_pilots(:,n))];
                An =  amp_symbol * U * D * diag(pilots(n,:)); %An
                A_bar = [real(An),-imag(An);imag(An),real(An)];   %��An��
                A_sum = A_sum + A_bar.' * A_bar;  %�˻���ǰ�벿��
                Ay = Ay + A_bar.' * yn; %�˻����벿��
            end
            theta_temp = A_sum \ Ay;   %ǰ�벿����˺�벿��
            theta(:,1) = theta_temp(1,1) + theta_temp(2,1) * 1i; %�˵�λ��͸�����λ��
            H_evd = U * diag(theta(:,1)); %ȡ�Խǲ���Ϊģ������
            
            
            %SVD
            [W,~,~] = svd(r_x);
            Ud = W(:,1);
            H_svd = Ud * Ud' * H_ls;



        
        %% MSE����
        csi_error = zeros(1,3);

        sum_e_H = zeros(1,3);
        for k = 1 : i_ant
            sum_e_H(1,1) = sum_e_H(1,1) + (H(k,1)-H_ls(k,1)) * (H(k,1)-H_ls(k,1))';
            sum_e_H(1,2) = sum_e_H(1,2) + (H(k,1)-H_evd(k,1)) * (H(k,1)-H_evd(k,1))';
            sum_e_H(1,3) = sum_e_H(1,3) + (H(k,1)-H_svd(k,1)) * (H(k,1)-H_svd(k,1))';
        end
        csi_error = csi_error + sum_e_H;

        csi_mse = csi_mse + csi_error / i_ant;

%         %% ������������
% 
%                 cova1 = zeros(i_ant,i_ant);
%                 R1 = zeros(i_ant,i_ant*Nt);
%                 Q1 = zeros(Nd,i_ant*Nt);
%                 Djj = diag(D);
%                 b = [pilots;symbol(:,(j-1)*Nd+1:j*Nd).'];
%                 R_e = (amp_use * b * D)';
%                     for k = 1 : i_ant
%                         R1((k-1)*K+1:k*K,(k-1)*Nt+1:k*Nt) = R_e;
%                         Q_e = (H(k,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K) * D(:,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K))';%amp_use not needed
%                         for l = 1 : Nd
%                             Q1((l-1)*K+1:l*K,(k-1)*Nt+tao+l) = Q_e;
%                         end
%                     end
%                     noise_va = 0;
%                     noise_va = noise_va + 1;
%                     Q_QH = Q1*Q1';
%                     R_RH = zeros(K*i_ant,K*i_ant);
%                     R_QH = zeros(K*i_ant,K*Nd);
%                     for k = 1 : i_ant
%                         R_QH((k-1)*K+1:k*K,:) = R1((k-1)*K+1:k*K,(k-1)*Nt+1:k*Nt) * (Q1(:,(k-1)*Nt+1:k*Nt))';
%                         R_RH((k-1)*K+1:k*K,(k-1)*K+1:k*K) = R1((k-1)*K+1:k*K,(k-1)*Nt+1:k*Nt) * (R1((k-1)*K+1:k*K,(k-1)*Nt+1:k*Nt))';
%                     end
%                     RPR1 = R_RH - R_QH / Q_QH * R_QH';
%                     cova1 = cova1 + noise_va * eye(K*i_ant) / RPR1;
%                
%                 bound1 = bound1 + sum(diag(cova1)) / L / K / i_ant;
    end
    csi_mse = csi_mse / Num;
    %bound1 = real(bound1) / Num;
    %mse_store(ii,:) = [csi_mse,bound1];
    mse_store(ii,:) = csi_mse;
    sprintf('%d',ii)
end



figure(1)
semilogy(SNR,mse_store(:,1),'b-.o','MarkerSize',10)
hold on
semilogy(SNR,mse_store(:,2),'-m+','MarkerSize',10)
semilogy(SNR,mse_store(:,3),'g-d','MarkerSize',10)
%semilogy(SNR,mse_store(:,4),'-r*')
grid on ;
legend('LS','EVD','SVD','Location','Northeast');
xlabel('SNR(dB)','Fontsize',10,'Fontname','Times')
ylabel('MSE','Fontsize',10,'Fontname','Times')
gtext('������Ŀ=128 ������Ŀ=100')




% figure(1)
% semilogy(SNR,mse_store(:,1),'b-.o','MarkerSize',10)
% % hold on
% % semilogy(SNR,mse_store(:,2),'-m+')
% % semilogy(SNR,mse_store(:,3),'g-d','MarkerSize',10)
% % semilogy(SNR,mse_store(:,4),'-r*')
% grid on ;
% legend('������Ŀ=128�����ݷ��ų���=100');
% xlabel('SNR(dB)','Fontsize',10,'Fontname','Times')
% ylabel('MSE','Fontsize',10,'Fontname','Times')
