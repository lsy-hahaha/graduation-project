clear;
close all;

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
shadow_amp = lognrnd(mu,sigma);
beta1 = shadow_amp  * ((rc*0.01).^(-0.5*gamma)); %用于求功率
amp_use = amp / beta1;

a = linspace(0.1,1,20);        % Interferance leval value.
beta111 = 1;

NUM_ITER = 1000;

% Generate pilot signals.
        %% 生成导频信号 
        S = zeros(tao,K);%Chu 
        if 0 == mod(tao,2)
            for k = 1 : K
                for l = 0 : tao-1
                    ll = mod(l+k-1, tao);
                    S(l+1,k) = exp(1i*pi*ll*ll/tao);
                end
            end
        else
            for k = 1 : K
                for l = 0 : tao-1
                    ll = mod(l+k-1, tao);
                    S(l+1,k) = exp(1i*pi*ll*(ll+1)/tao);
                end
            end
        end

%存储误差
ls_error_vec = zeros(1,length(a));
mmse_error_vec = zeros(1,length(a));

for a_idx=1:1:length(a)    
    ls_error = 0;
    mmse_error = 0;

    for numIter = 1:1:NUM_ITER
        
        beta_sum = 0;
        receive_pilots = zeros(i_ant,tao);
        Gil = zeros(i_ant,K,L);
        for l=1:1:L
            
            % Generate channels.
            beta = a(a_idx);
            if(l == 1)
                beta = 1;
            end
        
            %生成信道

            D = sqrt(beta)*eye(K);
            Gil(:,:,l) = (1/sqrt(2))*complex(randn(i_ant,K),randn(i_ant,K))*D;
            
            % Summation of all channels.
            receive_pilots = receive_pilots + Gil(:,:,l)*(S');
            
            % Summation of betas.
            beta_sum = beta_sum + beta;
            
        end
        
        % Factor.
        epsilon11 = (beta_sum + 1/(amp_use*tao));
        
        % Apply squared pilot power.
        receive_pilots = sqrt(amp_use)*receive_pilots;
        
        % Generate noise.
        W1 = (1/sqrt(2))*complex(randn(i_ant,tao),randn(i_ant,tao));
        
        % received pilot symbols at BS.
        Y1 = receive_pilots + W1;
        
        % ******* LS ********
        Z11_ls = (1/(sqrt(amp_use)*tao))*Y1*S(:,1);       
        ls_error = ls_error + ((Z11_ls'-Gil(:,1,1)')*(Z11_ls-Gil(:,1,1)));
        
        % ******* MMSE ********
        Z11_mmse = (beta111/epsilon11)*Z11_ls;
        mmse_error = mmse_error + ((Z11_mmse'-Gil(:,1,1)')*(Z11_mmse-Gil(:,1,1)));
        

    end
    
    fprintf('a: %d\n',a(a_idx));
    
    %LS
    ls_error = ls_error./(i_ant * NUM_ITER);
    ls_error_vec(a_idx) = ls_error;
    
    
    %MMSE
    mmse_error = mmse_error./(i_ant * NUM_ITER);
    mmse_error_vec(a_idx) = mmse_error;
    

end

fontSize = 10;

%fdee_figure = figure;
hold on;
semilogy(a,real(mmse_error_vec),'m-*','MarkerSize',8);
semilogy(a,real(ls_error_vec),'b-o','MarkerSize',8);

hold off
grid on;
xlabel('a')
ylabel('MSE')
legend('MMSE','LS', 'Location','northwest');

