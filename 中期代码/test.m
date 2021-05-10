clear;
close all;

%��С�����û�

K = 3;
Nd = 100;%���ų���

tao = 3;%��Ƶ����
r = 500;%С�����뾶
rc = r * 0.8;
rh = 100;%С����С�뾶
ra = rc / rh - 1;
gamma = 3.8;%·��˥��ָ��
mu = 0;%��Ӱ˥���ֵ
sigma = 10^0.8;%��Ӱ˥�䷽��
i_ant = 1000;%���������Ŀ
Num = 10;%����
Nt = tao + Nd;%�ܳ���

H = zeros(i_ant,1);
G = zeros(i_ant,1);

for ii = 1 : i_ant

        shadow_amp = lognrnd(mu,sigma);%������Ӱ˥��

       %�����ŵ�
        H = 1 / sqrt(2) * (randn(ii,K)+1i*randn(ii,K));
        D =  shadow_amp  * ((rc*0.01).^(-0.5*gamma)); %��߶�˥�� ��ֵ
        G = H * D;
        
        x(ii) = trace((H'*H))/ii/K;
        x1(ii) = abs(H(:,1)'*H(:,2)/ii);

end
plot(real(x));
hold on; 
grid on;
plot(real(x1));
xlabel('������Ŀ');
ylabel('�ڻ�');
legend('|Hi|^2/M','Hi��*Hj/M');
