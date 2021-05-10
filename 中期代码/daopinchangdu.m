clear;
close all;


%% ���� ֻ��LS��MMSE

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
Num = 10;%����
Nt = tao + Nd;%�ܳ���
SNR = 15;%�����
amp = 10 ^ (SNR*0.05) / sqrt(i_ant);
amp = 10 ^ (SNR*0.05) / sqrt(i_ant);

SNR_n = 7;
SNR = [0;5;10;15;20;25;30];
D = zeros(K,K*L*L);
H = zeros(i_ant,K*L*L);
G = zeros(i_ant,K*L*L);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mse_store = zeros(SNR_n, 2);

for ii = 1 : SNR_n
    amp = 10 ^ (SNR(ii,1)*0.05) / sqrt(i_ant);
    bound = 0;
    bound1 = 0;
    csi_mse = zeros(1,2);
   for jj = 1 : Num
       
        beta1 = a;
        if(j == 1)
           beta1 = 1;
        end
        betaMatrix = sqrt(beta1)*eye(K);



       %�����ŵ�
        for l1 = 1 : L%������վ
            for l2 = 1 : L
                H(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) = 1 / sqrt(2) * (randn(i_ant,K)+1i*randn(i_ant,K));
                D(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) = 
                G(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) = H(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) * D(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K);
            end
        end
        beta = shadow_amp  * ((rc*0.01).^(-0.5*gamma)); %��������
        amp_use = amp / beta;

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
        H_mmse = zeros(i_ant,K*L);
        H_ls = zeros(i_ant,K*L);
        
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
        end

            %MMSE
            %�ü򻯵����۷���
            beta_sum = zeros(1,1);
            for w = 1:K*L
                if w == 3 || 6 || 9
                    w_use =3;
                else
                w_use = mod(w,3);
                end
                
             beta_11 = D(w_use,(w_use-1)*L*K + (w_use-1)*K + w_use); 
             beta_sum = trace(D(:,3*w-2:3*w)) + 1/(tao*(amp)^2);
             beta_mmse = beta_11/beta_sum;
             H_mmse(:,w) = beta_mmse * H_ls(:,w);
            end
                
%                 if w == 1 || 2 || 3
%                     w1 = w;
%                     if w == 4 || 5 || 6
%                         w1 = w+6;
%                     if w == 7 || 8 || 9
%                         w1 = w+12;
%                     end
%                     end
%                 end
%                 
%                if w == 1 || 2 || 3
%                     w2 = w;
%                     if w == 4 || 5 || 6
%                         w1 = w+6;
%                     if w == 7 || 8 || 9
%                         w2 = w1+6;
%                     end
%                     end
%                 end
%                 
%                 beta_11 = D(w_use,w2); 
%                 beta_sum = D(w_use,w1) +D(w_use,w1+3) +D(w_use,w1+6) + 1/(tao*(amp_use)^2);
%                 beta_mmse = beta_11/beta_sum;
%                 H_mmse(:,w) = beta_mmse * H_ls(:,w);
%             end

%            %% MMSE ������ļ��㷽��
%             %��rh �ŵ����������ؾ���
%             for j1 = 1:L
%                 R_h = R_h + H_ls(:,(j1-1)*K+1:j1*K)'*H_ls(:,(j1-1)*K+1:j1*K);
%             end
%             H_ni = pilots' * R_h * pilots + eye(K) * 
 
            

        
        %% MSE����
        csi_error = zeros(1,2);
        for j = 1 : L
            sum_e_H = zeros(1,2);
            e_H = zeros(1,1);
            for k = 1 : i_ant
                sum_e_H(1,1) = sum_e_H(1,1) + (H(k,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K)-H_ls(k,(j-1)*K+1:j*K)) * (H(k,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K)-H_ls(k,(j-1)*K+1:j*K))';
                sum_e_H(1,2) = sum_e_H(1,2) + (H(k,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K)-H_mmse(k,(j-1)*K+1:j*K)) * (H(k,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K)-H_mmse(k,(j-1)*K+1:j*K))';
%                 e_H = e_H + (H(k,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K))* (H(k,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K))';
%                 sum_e_H  = sum_e_H./e_H;
            end
            csi_error = csi_error + sum_e_H;
        end
        csi_mse = csi_mse + csi_error / L / K / i_ant;


    end
    csi_mse = csi_mse / Num;
    mse_store(ii,:) = [csi_mse];
    sprintf('%d',ii)
end



figure(1)
semilogy(SNR,mse_store(:,1),'b-.o','MarkerSize',8)
hold on
semilogy(SNR,mse_store(:,2),'-m+','MarkerSize',8)
grid on ;
legend('LS','MMSE','Location','Northeast');
xlabel('SNR(dB)','Fontsize',10,'Fontname','Times')
ylabel('MSE','Fontsize',10,'Fontname','Times')
gtext('������Ŀ=128 ������Ŀ=100')

