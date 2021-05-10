clear;
close all;

%% ������С������ŵ����

L = 3;%С����Ŀ
Nd = 100;%���ų���
K = 3;%�û���Ŀ
tao = K;%��Ƶ����
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
SNR = 15;%�����
amp = 10 ^ (SNR*0.05) / sqrt(i_ant);
%%ÿ����վ��λ��
% base(1:7,1) = [0;(1i * 2 * rc);(sqrt(3) * rc + 1i * rc);(sqrt(3) * rc - 1i * rc);(-1i * 2 * rc);(-sqrt(3) * rc - 1i * rc);(-sqrt(3) * rc + 1i * rc);];
% beta = 1;
SNR_n = 11;
SNR = [0;2;4;6;8;10;12;14;16;18;20];
D = zeros(K,K*L*L);
H = zeros(i_ant,K*L*L);
G = zeros(i_ant,K*L*L);
a = 0.05;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mse_store = zeros(SNR_n, 5);

for ii = 1 : SNR_n
    amp = (10 ^ (SNR(ii,1)./10)) / (tao*(1 + a*(L-1)));
    amp_use = sqrt(amp);
    bound = 0;
    bound1 = 0;
    csi_mse = zeros(1,5);
   for jj = 1 : Num
        %shadow_amp = lognrnd(mu,sigma);%������Ӱ˥��
        %%�����û��ն˾��ȷֲ�
%         dis(1:K,1:3) = (rem(rand(K,3) * ra, ra) + 1) * rh;
%         ang(1:K,1:3) = rand(K,3) * 2 * pi;
%         pos(1:K,1:3) = dis .* (exp(1i * ang));
%         pos(:,2) = pos(:,2) + base(2,1);
%         pos(:,3) = pos(:,3) + base(3,1);

       %�����ŵ�
       beta_sum = 0;
       for j = 1:L
       beta1 = a;
        if(j == 1)
           beta1 = 1;
        end
        betaMatrix = sqrt(beta1)*eye(K);
        
        beta111 = 1;
        beta_sum = beta_sum + beta1;
        epsilon11 = beta_sum + 1/(amp*tao);
       end
       
%         for l1 = 1 : L%������վ
%             for l2 = 1 : L
%                 H(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) = 1 / sqrt(2) * (randn(i_ant,K)+1i*randn(i_ant,K));
%                 D(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) = shadow_amp  * diag(((abs(pos(:,l2)-base(l1,1))*0.01).^(-0.5*gamma)));
%                 G(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) = H(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) * D(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K);
%             end
%         end

        for l1 = 1 : L%������վ
            for l2 = 1 : L
                H(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) = 1 / sqrt(2) * (randn(i_ant,K)+1i*randn(i_ant,K));
                D(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) =  betaMatrix;
                G(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) = H(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) * D(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K);
            end
        end
%         beta = shadow_amp  * ((rc*0.01).^(-0.5*gamma)); %��������
%         amp_use = amp / beta;

        %% ���ɵ��Ʒ��ţ�ÿN�ж�Ӧһ��С���еķ���
        symbol = sign(randn(K,L*Nd));%BPSK����
        for k = 1 : K
            for l = 1 : L*Nd
                if 0 == symbol(k,l)
                    symbol(k,l) = 1;%ֻ�С�1
                else
                end
            end
        end
        
       %% ���ɽ����ź� ���������С��֮��ĸ��ŵĻ� ֻȡG��1:3 13:15��25:27
        %һ��ʼ���ǵ��Ǹ��źʹ����źŹ�������ͬ��
        receive_symbol = zeros(i_ant,Nd*L); %�洢�����ź�
        noise =  (randn(i_ant,Nd*L)+1i*randn(i_ant,Nd*L)) / sqrt(2);%��������
        receive_symbol = receive_symbol + noise; %����
        for j = 1 : L
           for l = 1 : L
               if j == l 
                Gjl = G(:,(j-1)*L*K+(l-1)*K+1:(j-1)*L*K+l*K);
               else
                   Gjl = zeros(i_ant,L);
               end
                receive_symbol(:,(j-1)*Nd+1:j*Nd) = receive_symbol(:,(j-1)*Nd+1:j*Nd) + amp_use * Gjl * symbol(:,(j-1)*Nd+1:j*Nd);
           end
        end
        

        %% ���ɵ�Ƶ�ź� 
        pilots = zeros(tao,K);%Chu 
        if 0 == mod(tao,2)
            for k = 1 : K
                for l = 0 : tao-1
                    ll = mod(l+k-1, tao);
                    pilots(l+1,k) = exp(1i*pi*ll*ll/tao);
                end
            end
        else
            for k = 1 : K
                for l = 0 : tao-1
                    ll = mod(l+k-1, tao);
                    pilots(l+1,k) = exp(1i*pi*ll*(ll+1)/tao);
                end
            end
        end
        R_p = pilots.' * conj(pilots);
        %% ���յ�Ƶ�ź�
        receive_pilots = zeros(i_ant,tao*L);
        noise1 =  (randn(i_ant,tao*L)+1i*randn(i_ant,tao*L)) / sqrt(2);
        receive_pilots = receive_pilots + noise1;
        for j = 1 : L
            for l = 1 : L
                if j == l
                Gjl = G(:,(j-1)*L*K+(l-1)*K+1:(j-1)*L*K+l*K);
                else
                    Gjl = zeros(i_ant,L);
                end
                receive_pilots(:,(j-1)*tao+1:j*tao) = receive_pilots(:,(j-1)*tao+1:j*tao) + amp_use * Gjl  * pilots.';
            end
        end
        
        
        %% �ŵ�����
        H_svd = zeros(i_ant,K*L);
        H_ls = zeros(i_ant,K*L);
        theta = zeros(K,L);
        H_evd = zeros(i_ant,K*L);
        U = zeros(i_ant,K);
        
        for j = 1 : L
            r_x = zeros(i_ant,i_ant);
            for k = 1 : tao
                r_x = r_x + receive_pilots(:,(j-1)*tao+k) * receive_pilots(:,(j-1)*tao+k)';
            end
            for k = 1 : Nd
                r_x = r_x + receive_symbol(:,(j-1)*Nd+k) * receive_symbol(:,(j-1)*Nd+k)';
            end
            r_x = r_x / Nt;   %������źŵ�Э������� ��Ƶ������

            %LS
            Hpi = receive_pilots(:,(j-1)*tao+1:j*tao) * conj(pilots) / R_p / amp_use;
            H_ls(:,(j-1)*K+1:j*K) = Hpi / D(:,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K);
            
            %MMSE
            H_mmse(:,(j-1)*K+1:j*K) = (beta111/epsilon11)*H_ls(:,(j-1)*K+1:j*K);

            
            %EVD 
            [U_all,~] = eig(r_x);
            U = U_all(:,1:K);
            A_sum = zeros(2*K,2*K);
            Ay = zeros(2*K,1);
            for n = 1 : tao
                yn = [real(receive_pilots(:,(j-1)*tao+n));imag(receive_pilots(:,(j-1)*tao+n))];
                An =  amp_use * U * (D(:,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K)) * diag(pilots(:,n)); %An
                A_bar = [real(An),-imag(An);imag(An),real(An)];   %��An��
                A_sum = A_sum + A_bar.' * A_bar;  %�˻���ǰ�벿��
                Ay = Ay + A_bar.' * yn; %�˻����벿��
            end
            theta_temp = A_sum \ Ay;   %ǰ�벿����˺�벿��
            theta(:,j) = theta_temp(1:K,1) + theta_temp(K+1:2*K,1) * 1i; %�˵�λ��͸�����λ��
            H_evd(:,(j-1)*K+1:j*K) = U * diag(theta(:,j)); %ȡ�Խǲ���Ϊģ������
            
            
            %SVD
            [W,~,~] = svd(r_x);
            Ud = W(:,1:K);
            H_svd(:,(j-1)*K+1:j*K) = Ud * Ud' * H_ls(:,(j-1)*K+1:j*K);
            
            % ����8�еķ���
            Uc = W(:,K+1:i_ant);
            R_u = zeros(i_ant,i_ant);
            for l = 1 : i_ant-K
                R_u = R_u + Uc(:,l) * Uc(:,l)';
            end
            for k = 1 : K
                sum_H_e_pilot = zeros(i_ant,1);
                R_temp = R_u;
                for j1 = 1 : K
                    R_temp = R_temp +  (H_ls(:,(j-1)*K+j1) * H_ls(:,(j-1)*K+j1)');
                end
                H_8(:,(j-1)*K+k) =  i_ant * (R_temp \ H_ls(:,(j-1)*K+k));
            end
            
            
            % SVD+ILSP
            


        end
        
        
        % SVD+ILSP
        
        
        
        %% MSE����
        csi_error = zeros(1,5);
        for j = 1 : L
            sum_e_H = zeros(1,5);
            e_H = zeros(1,1);
            for k = 1 : i_ant
                sum_e_H(1,1) = sum_e_H(1,1) + (H(k,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K)-H_ls(k,(j-1)*K+1:j*K)) * (H(k,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K)-H_ls(k,(j-1)*K+1:j*K))';
                sum_e_H(1,2) = sum_e_H(1,2) + (H(k,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K)-H_evd(k,(j-1)*K+1:j*K)) * (H(k,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K)-H_evd(k,(j-1)*K+1:j*K))';
                sum_e_H(1,3) = sum_e_H(1,3) + (H(k,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K)-H_svd(k,(j-1)*K+1:j*K)) * (H(k,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K)-H_svd(k,(j-1)*K+1:j*K))';
                sum_e_H(1,4) = sum_e_H(1,4) + (H(k,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K)-H_8(k,(j-1)*K+1:j*K)) * (H(k,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K)-H_8(k,(j-1)*K+1:j*K))';
                sum_e_H(1,5) = sum_e_H(1,5) + (H(k,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K)-H_mmse(k,(j-1)*K+1:j*K)) * (H(k,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K)-H_mmse(k,(j-1)*K+1:j*K))';
                e_H = e_H + (H(k,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K))* (H(k,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K))';
                sum_e_H  = sum_e_H./e_H;
            end
            csi_error = csi_error + sum_e_H;
        end
        %csi_mse = csi_mse + csi_error / L / K / i_ant;
        csi_mse = csi_mse + csi_error;

        %% ������������

%                 cova1 = zeros(K*i_ant,K*i_ant);
%                 R1 = zeros(K*i_ant,i_ant*Nt);
%                 Q1 = zeros(K*Nd,i_ant*Nt);
%                 for j = 1 : L
%                     Djj = diag(D(:,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K));
%                     b = [pilots;symbol(:,(j-1)*Nd+1:j*Nd).'];
%                     R_e = (amp_use * b * D(:,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K))';
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
%                 end
%                 for i = 1:K
%                     for j = 1:L
%                         e_H = e_H + (H(k,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K))* (H(k,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K))';
%                     end
%                 end
%                 bound1 = bound1 + sum(diag(cova1)) / e_H;
    end
%     csi_mse = csi_mse / Num;
%     %bound = real(bound) / Num;
%     bound1 = real(bound1) / Num;
%     mse_store(ii,:) = [csi_mse,bound1];
    
    csi_mse = csi_mse / Num;
    mse_store(ii,:) = [csi_mse];
    sprintf('%d',ii)
end



figure(1)
semilogy(SNR,mse_store(:,1),'b-.o','MarkerSize',7)
hold on
semilogy(SNR,mse_store(:,2),'-m+','MarkerSize',7)
semilogy(SNR,mse_store(:,3),'g-d','MarkerSize',7)
semilogy(SNR,mse_store(:,4),'k-^','MarkerSize',7)
semilogy(SNR,mse_store(:,5),'-r*','MarkerSize',7)
grid on ;
legend('LS','EVD','SVD','����8','MMSE','Location','Northeast');
xlabel('SNR(dB)','Fontsize',10,'Fontname','Times')
ylabel('MSE','Fontsize',10,'Fontname','Times')
gtext('������Ŀ=128 ������Ŀ=100')

