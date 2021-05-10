clear;
close all;

%% 不考虑小区间干扰的情况

L = 3;%小区数目
Nd = 100;%符号长度
K = 3;%用户数目
tao = K;%导频长度
r = 500;%小区最大半径
rc = r * 0.8;
rh = 100;%小区最小半径
ra = rc / rh - 1;
gamma = 3.8;%路径衰减指数
mu = 0;%阴影衰落均值
sigma = 10^0.8;%阴影衰落方差
i_ant = 128;%最大天线数目
Num = 100;%迭代
Nt = tao + Nd;%总长度
SNR = 15;%信噪比
amp = 10 ^ (SNR*0.05) / sqrt(i_ant);
ant_n = 7;
ant_s = [40;80;120;160;200;240;280];%天线数目范围

a = 0.05;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mse_store = zeros(ant_n, 5);

for ii = 1 : ant_n
    i_ant = ant_s(ii,1);
    
    D = zeros(K,K*L*L);
    H = zeros(i_ant,K*L*L);
    G = zeros(i_ant,K*L*L);
    
    amp = (10 ^ (SNR./10)) / (tao*(1 + a*(L-1)));

    amp_use = sqrt(amp);
    bound = 0;
    bound1 = 0;
    csi_mse = zeros(1,5);
   for jj = 1 : Num
        %shadow_amp = lognrnd(mu,sigma);%定义阴影衰落
        %%假设用户终端均匀分布
%         dis(1:K,1:3) = (rem(rand(K,3) * ra, ra) + 1) * rh;
%         ang(1:K,1:3) = rand(K,3) * 2 * pi;
%         pos(1:K,1:3) = dis .* (exp(1i * ang));
%         pos(:,2) = pos(:,2) + base(2,1);
%         pos(:,3) = pos(:,3) + base(3,1);

       %生成信道
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
       
%         for l1 = 1 : L%遍历基站
%             for l2 = 1 : L
%                 H(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) = 1 / sqrt(2) * (randn(i_ant,K)+1i*randn(i_ant,K));
%                 D(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) = shadow_amp  * diag(((abs(pos(:,l2)-base(l1,1))*0.01).^(-0.5*gamma)));
%                 G(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) = H(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) * D(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K);
%             end
%         end

        for l1 = 1 : L%遍历基站
            for l2 = 1 : L
                H(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) = 1 / sqrt(2) * (randn(i_ant,K)+1i*randn(i_ant,K));
                D(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) =  betaMatrix;
                G(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) = H(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) * D(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K);
            end
        end
%         beta = shadow_amp  * ((rc*0.01).^(-0.5*gamma)); %用于求功率
%         amp_use = amp / beta;

        %% 生成调制符号，每N列对应一个小区中的符号
        symbol = sign(randn(K,L*Nd));%BPSK调制
        for k = 1 : K
            for l = 1 : L*Nd
                if 0 == symbol(k,l)
                    symbol(k,l) = 1;%只有±1
                else
                end
            end
        end
        
       %% 生成接收信号 如果不考虑小区之间的干扰的话 只取G的1:3 13:15和25:27
        %一开始考虑的是干扰和传输信号功率是相同的
        receive_symbol = zeros(i_ant,Nd*L); %存储接收信号
        noise =  (randn(i_ant,Nd*L)+1i*randn(i_ant,Nd*L)) / sqrt(2);%生成噪声
        receive_symbol = receive_symbol + noise; %加噪
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
        

        %% 生成导频信号 
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
        %% 接收导频信号
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
        
        
        %% 信道估计
        H_svd = zeros(i_ant,K*L);
        H_ls = zeros(i_ant,K*L);
        H_mmse = zeros(i_ant,K*L);
        H_8 = zeros(i_ant,K*L);
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
            r_x = r_x / Nt;   %求接收信号的协方差矩阵 导频加数据

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
                A_bar = [real(An),-imag(An);imag(An),real(An)];   %求An拔
                A_sum = A_sum + A_bar.' * A_bar;  %乘积项前半部分
                Ay = Ay + A_bar.' * yn; %乘积项后半部分
            end
            theta_temp = A_sum \ Ay;   %前半部分逆乘后半部分
            theta(:,j) = theta_temp(1:K,1) + theta_temp(K+1:2*K,1) * 1i; %乘单位阵和复数单位阵
            H_evd(:,(j-1)*K+1:j*K) = U * diag(theta(:,j)); %取对角部分为模糊矩阵
            
            
            %SVD
            [W,~,~] = svd(r_x);
            Ud = W(:,1:K);
            H_svd(:,(j-1)*K+1:j*K) = Ud * Ud' * H_ls(:,(j-1)*K+1:j*K);
            
            % 文献8中的方法
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
        
        
        
        %% MSE计算
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
                %e_H = e_H + (H(k,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K))* (H(k,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K))';
                %sum_e_H  = sum_e_H./e_H;
            end
            csi_error = csi_error + sum_e_H;
        end
        csi_mse = csi_mse + csi_error / L / K / i_ant;
        %csi_mse = csi_mse + csi_error;

        %% 克拉美罗下限

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
semilogy(ant_s,mse_store(:,1),'b-.o','MarkerSize',7)
hold on
%semilogy(SNR,mse_store(:,2),'-m+','MarkerSize',7)
%semilogy(SNR,mse_store(:,3),'g-d','MarkerSize',7)
%semilogy(SNR,mse_store(:,4),'k-^','MarkerSize',7)
semilogy(ant_s,mse_store(:,5),'-r*','MarkerSize',7)
grid on ;
legend('LS','MMSE','Location','Northeast');
xlabel('number of antennas','Fontsize',10,'Fontname','Times')
ylabel('MSE','Fontsize',10,'Fontname','Times')
axis([40 280 0.1 1])
gtext('SNR=15dB 符号数目=100')

