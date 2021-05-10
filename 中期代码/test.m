clear;
close all;

%单小区单用户

K = 3;
Nd = 100;%符号长度

tao = 3;%导频长度
r = 500;%小区最大半径
rc = r * 0.8;
rh = 100;%小区最小半径
ra = rc / rh - 1;
gamma = 3.8;%路径衰减指数
mu = 0;%阴影衰落均值
sigma = 10^0.8;%阴影衰落方差
i_ant = 1000;%最大天线数目
Num = 10;%迭代
Nt = tao + Nd;%总长度

H = zeros(i_ant,1);
G = zeros(i_ant,1);

for ii = 1 : i_ant

        shadow_amp = lognrnd(mu,sigma);%定义阴影衰落

       %生成信道
        H = 1 / sqrt(2) * (randn(ii,K)+1i*randn(ii,K));
        D =  shadow_amp  * ((rc*0.01).^(-0.5*gamma)); %大尺度衰落 定值
        G = H * D;
        
        x(ii) = trace((H'*H))/ii/K;
        x1(ii) = abs(H(:,1)'*H(:,2)/ii);

end
plot(real(x));
hold on; 
grid on;
plot(real(x1));
xlabel('天线数目');
ylabel('内积');
legend('|Hi|^2/M','Hi＇*Hj/M');
