% 论文仿真原型程序,仿真五种（MRC,ZF,ZF-SIC,MMSE,MMSE-SIC） Vblast接收的检测性能，绘制误比特率～信噪比曲线。 
% 发端初始化================================================================ 
% 发射天线数tx,接收天线数rx,发射矩阵长度L(帧长)
clear all;
clc;

tx=4;rx=4;L=10000; 
Modulation='QPSK'; 
EbN0=[0:5:20]; 
B=30000;Ts=1/24300; 
% 建立EbN0与SNR之间的换算关系 
SNR=EbN0-10*log10(Ts*B); 
% 信源A 
A=randi([0,3],tx*L,1); 
% 经过QPSK调制的V-Blast发射矩阵X 
X1=reshape(A,tx,L);
X=zeros(tx,L);
for k=1:L
    X(:,k)=pskmod(X1(:,k),4);
end

% 信道传输=================================================================

% 快衰落Rayleigh信道H 
H=sqrt(1/2)*(randn(rx,tx,L)+i*randn(rx,tx,L)); 
% 均值为0方差为1的复高斯白噪声n 
n=sqrt(1/2)*(randn(rx,L)+i*randn(rx,L)); 
% 未叠加噪声的接收信号R 
R=zeros(rx,L); 
for k=1:L 
    R(:,k)=sqrt(1/tx)*H(:,:,k)*X(:,k); 
end 

% 检测 

%MRC=======================================================================
disp('MRC')
berm=[];

for m=SNR
    m
    snr=10^(m/10);
    R_noised=R+sqrt(1/snr)*n; 
    x=[];
    for t=1:L
        r=R_noised(:,t);
        HH=H(:,:,t);
        G=HH';
        y=G*r;
        xtemp=zeros(tx,1);
        for k=1:tx
            if imag(y(k))<real(y(k)) && imag(y(k))>-real(y(k))
                xtemp(k)=1;
            elseif imag(y(k))>real(y(k)) && imag(y(k))>-real(y(k))
                xtemp(k)=i;
            elseif imag(y(k))>real(y(k)) && imag(y(k))<-real(y(k))
                xtemp(k)=-1;
            elseif imag(y(k))<real(y(k)) && imag(y(k))<-real(y(k))
                xtemp(k)=-i;
            end
        end
        x=[x,xtemp];
    end
    % 从x求A的估计a 
    x1=zeros(tx,L);
    for k=1:L
        x1(:,k)=pskdemod(x(:,k),4);
    end
    a=reshape(x1,tx*L,1);
    % 比较A和a计算错值率temp_ber 
    [errbit,temp_ber]=biterr(A,a,2); 
    berm=[berm,temp_ber];
end
figure
semilogy(EbN0,berm,'*- g')
hold on
%ZF========================================================================
disp('ZF');
berz=[];
% 在不同的信噪比下计算ZF接收机误比特率berz 
for m=SNR 
    m
    % 每个子信道的平均信噪比为snr的接收信号R_noised 
    snr=10^(m/10);
    R_noised=R+sqrt(1/snr)*n; 
    x=[];
    % 逐时隙对接收符号矢量进行检测，合并得到一帧发射矩阵X的估计x 
    for t=1:L
        r=R_noised(:,t);
        % 迫零矩阵G 
        G=pinv(H(:,:,t));
        y=G*r;
        % QPSK判决
        xtemp=zeros(tx,1);
        for k=1:tx
            if imag(y(k))<real(y(k)) && imag(y(k))>-real(y(k))
                xtemp(k)=1;
            elseif imag(y(k))>real(y(k)) && imag(y(k))>-real(y(k))
                xtemp(k)=i;
            elseif imag(y(k))>real(y(k)) && imag(y(k))<-real(y(k))
                xtemp(k)=-1;
            elseif imag(y(k))<real(y(k)) && imag(y(k))<-real(y(k))
                xtemp(k)=-i;
            end
        end
        x=[x,xtemp];
    end
    % 从x求A的估计a 
    x1=zeros(tx,L);
    for k=1:L
        x1(:,k)=pskdemod(x(:,k),4);
    end
    a=reshape(x1,tx*L,1);
    % 比较A和a计算错值率temp_ber 
    [errbit,temp_ber]=biterr(A,a,2); 
    berz=[berz,temp_ber];
end
semilogy(EbN0,berz,'o- b')
% ZF-SIC(ordered)==========================================================
disp('ZF-SIC')
berzs=[];

for m=SNR
    m
    snr=10^(m/10);
    R_noised=R+sqrt(1/snr)*n; 
    x=[];
    for t=1:L
        r=R_noised(:,t);
        HH=H(:,:,t);
        G=pinv(HH);
        S=[1:tx];% S表示一个时隙内还未检测的符号的序号的集合
        xtemp=zeros(tx,1);
        % 逐发射天线进行检测 
        for k=1:tx
            % G范数最小的行是wki，它是G的第ki行 
            [wki,ki]=minnorm(G);
            % 判决统计量y 
            y=wki*r; 
            % QPSK判决
            if imag(y)<real(y) && imag(y)>-real(y)
                xtemp(S(ki))=1;
            elseif imag(y)>real(y) && imag(y)>-real(y)
                xtemp(S(ki))=i;
            elseif imag(y)>real(y) && imag(y)<-real(y)
                xtemp(S(ki))=-1;
            elseif imag(y)<real(y) && imag(y)<-real(y)
                xtemp(S(ki))=-i;
            end
            % SIC串行干扰抵消 
            r=r-sqrt(1/tx)*xtemp(S(ki))*H(:,S(ki),t);
            % 将已经检测的信号对应的信道矩阵的列删去 
            HH(:,ki)=[];
            % 已经检测过的信号的序号删去
            S(ki)=[];
            G=pinv(HH); 
        end
        x=[x,xtemp];
    end
    % 从x求A的估计a 
    x1=zeros(tx,L);
    for k=1:L
        x1(:,k)=pskdemod(x(:,k),4);
    end
    a=reshape(x1,tx*L,1);
    % 比较A和a计算错值率temp_ber 
    [errbit,temp_ber]=biterr(A,a,2); 
    berzs=[berzs,temp_ber];
end
semilogy(EbN0,berzs,'o- r')

% MMSE====================================================================
disp('MMSE')
bermm=[];
 
for m=SNR
    m
    snr=10^(m/10);
    R_noised=R+sqrt(1/snr)*n;
    x=[];
    for t=1:L
        r=R_noised(:,t);
        HH=H(:,:,t);
        G=inv(HH'*HH+(1/snr)*eye(tx))*HH';
        y=G*r;
        xtemp=zeros(tx,1);
        for k=1:tx
            if imag(y(k))<real(y(k)) && imag(y(k))>-real(y(k))
                xtemp(k)=1;
            elseif imag(y(k))>real(y(k)) && imag(y(k))>-real(y(k))
                xtemp(k)=i;
            elseif imag(y(k))>real(y(k)) && imag(y(k))<-real(y(k))
                xtemp(k)=-1;
            elseif imag(y(k))<real(y(k)) && imag(y(k))<-real(y(k))
                xtemp(k)=-i;
            end
        end
        x=[x,xtemp];
    end
    % 从x求A的估计a 
    x1=zeros(tx,L);
    for k=1:L
        x1(:,k)=pskdemod(x(:,k),4);
    end
    a=reshape(x1,tx*L,1);
    % 比较A和a计算错值率temp_ber 
    [errbit,temp_ber]=biterr(A,a,2); 
    bermm=[bermm,temp_ber];
end
semilogy(EbN0,bermm,'s- b') 
% MMSE-SIC(ordered)========================================================
disp('MMSE-SIC')
bermms=[];

for m=SNR
    m
    snr=10^(m/10);
    R_noised=R+sqrt(1/snr)*n;   
    x=[];
    for t=1:L
        r=R_noised(:,t);
        HH=H(:,:,t); 
        G=inv(HH'*HH+(1/snr)*eye(tx))*HH';
        S=[1:tx];% S表示一个时隙内还未检测的符号的序号的集合
        xtemp=zeros(tx,1);
        for k=1:tx
            [wki,ki]=minnorm(G);
            % 判决统计量y 
            y=wki*r;
            % QPSK判决
            if imag(y)<real(y) && imag(y)>-real(y)
                xtemp(S(ki))=1;
            elseif imag(y)>real(y) && imag(y)>-real(y)
                xtemp(S(ki))=i;
            elseif imag(y)>real(y) && imag(y)<-real(y)
                xtemp(S(ki))=-1;
            elseif imag(y)<real(y) && imag(y)<-real(y)
                xtemp(S(ki))=-i;
            end
            % SIC串行干扰抵消 
            r=r-sqrt(1/tx)*xtemp(S(ki))*H(:,S(ki),t);
            % 将已经检测的信号对应的信道矩阵的列删去 
            HH(:,ki)=[];
            % 已经检测过的信号的序号删去
            S(ki)=[];
            G=inv(HH'*HH+(1/snr)*eye(tx-k))*HH';
        end
        x=[x,xtemp];
    end
    % 从x求A的估计a 
    x1=zeros(tx,L);
    for k=1:L
        x1(:,k)=pskdemod(x(:,k),4);
    end
    a=reshape(x1,tx*L,1);
    % 比较A和a计算错值率temp_ber 
    [errbit,temp_ber]=biterr(A,a,2); 
    bermms=[bermms,temp_ber];
end
semilogy(EbN0,bermms,'s- r')
grid on

legend('MRC','ZF','ZF-SIC','MMSE','MMSE-SIC'); 
xlabel('EbN0(dB)'); 
ylabel('误比特率'); 
title('blast检测比较')
        

