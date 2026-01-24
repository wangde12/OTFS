%Simultaneous
%Run 'main.m', and you will get the results shown in Fig. 8 and Fig. 9. ;
%'168R.mat' is proposed Zak sequence in Zak domain (DD domain) ;
% Contact email： pengxp@ysu.edu.cn
%Thanks to Raviteja Patchava et al. for sharing the procedures ('OTFS_channel_gen.m','OTFS_channel_output.m' and 'OTFS_modulation.m')
%They can be downloaded from *OTFS slides* : https://www.ecse.monash.edu.au/staff/eviterbo/OTFS-VTC18/index.html



clear
clc
close all
rng("shuffle")


N = 16;
% number of subcarriers
M = 8;
numSym = N*M;
Nfft=N*M;


load('168R.mat');
frame_num=3;%Three frames, the middle frame is the synchronization frame.

X1=zeros(16,8);
X1(1,1)=exp(0*2*pi*i/16);
X1(2,2)=exp(2*2*pi*i/16);
X1(3,7)=exp(12*2*pi*i/16);
X1(4,3)=exp(4*2*pi*i/16);
X1(5,6)=exp(10*2*pi*i/16);
X1(6,4)=exp(6*2*pi*i/16);
X1(7,8)=exp(14*2*pi*i/16);
X1(8,5)=exp(8*2*pi*i/16);
X1(9,1)=exp(4*2*pi*i/16);
X1(10,2)=exp(6*2*pi*i/16);
X1(11,7)=exp(0*2*pi*i/16);
X1(12,3)=exp(8*2*pi*i/16);
X1(13,6)=exp(4*2*pi*i/16);
X1(14,4)=exp(10*2*pi*i/16);
X1(15,8)=exp(2*2*pi*i/16);
X1(16,5)=exp(12*2*pi*i/16);
X1=X1/abs(X1(1,1));
po=sum(sum(abs(X1).^2));

% size of constellation
M_mod = 4;
M_bits = log2(M_mod);

twiddle_factor = exp(-1j*2*pi*(0:M-1).'*(0:N-1)/(N*M));  % 旋转因子

% 发射脉冲 / 接收脉冲
% Delta-like
g_tx = zeros(M,N);
g_tx(1,1) = 1;  % MATLAB索引从1开始
g_rx = conj(g_tx);  % 匹配滤波

% 构造 IDZT 和 DZT 的矩阵部分
F_inv = ifft(eye(N)) * sqrt(N);  % N x N, 对每列执行 IFFT*sqrt(N)
F = fft(eye(N)) / sqrt(N);       % N x N
I_F_inv=kron(eye(M), F_inv);
I_F=kron(eye(M), F);

% average energy per data symbol
eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));
X1=	sqrt(((M*N)*2)/po)*X1;
N_syms_perfram = N*M;
% number of bits per frame
N_bits_perfram = N*M*M_bits;

Fn1=dftmtx(N);
Fn=Fn1/norm(Fn1);

% st1 = OTFS_modulation(N,M,X1);%OTFS modulation

vec_X1=reshape(abs(X1.'),N*M,1);
[index_value,~]=find(vec_X1>0);
index_zero=setdiff(1:N*M,index_value);

R1=sum(sum(X1.*conj(X1)));  % R1=array_M*array_N:二维自相关函数在零移位时的值(自相关函数峰值位置对应的值)

X_td = (X1.').*conj(twiddle_factor).*(twiddle_factor);
X_tw = twisted_conv(X_td, g_tx, N, M);    % 发射端：扭曲卷积 + IDZT
% s = OTFS_modulation(N,M,(X));
s_mat=ifft(X_tw, [],2)*sqrt(N);
st1 = s_mat(:);

% KK2_=pccf(st1,st1);
% plot(abs(KK2_))

% KP=[4,1;4,2;8,1;8,2;12,1;12,2;4,1;4,2;8,1;8,2;12,1;12,2];% first 6 row is Zak sequence in this paper, later 6 row is random sequence.
KP=[6,1;6,1];
% KP=[4,1;4,1];
taps = KP(1,1);

% if value==1
%     KP=[4,1;4,1];
% else
%     KP=[4,2;4,2];
% end
% end
BL=1e2;%Monte Carlo num
SNR_dB =-10:2:20;
% SNR_dB =20;
cpl=32;%CP length
TO=zeros(length(SNR_dB),1);
R_all=[];
v_max=200; % km/h

ber_snr_1=zeros(BL,length(KP));
for cc=1:length(KP)
    for num1=1:length(SNR_dB)
        KK=0;
        KK_DD=0;
        KK_equalization_DD=0;
        amp=sqrt(2);
        sigma2=10^(-SNR_dB(num1)/10)*amp^2;

        % Ber_est  = 3; %
        % FramNum  = 100+10*10.^(Ber_est*num1/length(SNR_dB));
        % BL=floor(FramNum);

        ber_fram_true=zeros(BL,1);

        is_syn=zeros(BL,1);
        pos_syn=[];
     %   kmk=[];
       % for num2=1:BL
       ber_snr_equalization_DD=zeros(BL,1);
       ber_snr_DD=zeros(BL,1);
       ber_snr=zeros(BL,1);
       parfor num2=1:BL
            R_all=[];
            [delay_taps,Doppler_taps,chan_coef,f_max] = OTFS_channel_gen6(N,M,taps,v_max);      %% OTFS channel generation%%%%
            carrier_space=15*10^3;
            k_max=f_max*N/carrier_space;
            % DD信道矩阵
            H_eff = zeros(numSym,numSym);
            for k = 1:numSym
                e = zeros(numSym,1); e(k)=1;
                % 发射端卷积 + IDZT
                E_td=reshape(e,M,N).*twiddle_factor;
                E_tw=twisted_conv(E_td, g_tx,N,M);
                Y_E=ifft(E_tw, [],2)*sqrt(N);
                s_e = reshape(Y_E,[],1);
                r_e = OTFS_channel_output0(N,M,taps,delay_taps,Doppler_taps,chan_coef,0,s_e);
                % 接收端匹配卷积
                r_e_blocks = reshape(r_e,M,N);
                R_e = fft(r_e_blocks,[],2)/sqrt(N);
                R_e_matched = twisted_conv(R_e, conj(g_rx), N, M);
                R_e_matched = (R_e_matched) .* conj(twiddle_factor);
                H_eff(:,k) = R_e_matched(:);
            end
            HDD=H_eff;
            Wmmse_tradition=(HDD'*HDD+sigma2*eye(N*M))\HDD';

            Y=zeros(taps,3*(N*M+cpl)+2*cpl+2*cpl);
            for num3=1:3%frame_num
                data_info_bit = randi([0,1],N_bits_perfram,1);
                data_temp = bi2de(reshape(data_info_bit,N_syms_perfram,M_bits));
                x = qammod(data_temp,M_mod,'gray');
                X = reshape(x,M,N);
                if num3==2
                    X=X1.'.*conj(twiddle_factor);
                    if  cc>=2
                        data_info_bit = randi([0,1],N_bits_perfram,1);
                        data_temp = bi2de(reshape(data_info_bit,N_syms_perfram,M_bits));
                        x = qammod(data_temp,M_mod,'gray');
                        X = reshape(x,M,N);

                        % X_td = X.*(twiddle_factor);
                        % X_tw = twisted_conv(X_td.', g_tx, N, M);    % 发射端：扭曲卷积 + IDZT
                        % X_tw=X_tw.';
                        % s_mat=X_tw*Fn';
                        % st2 = s_mat(:);

                        % s_mat=X*Fn';
                        % st2 = s_mat(:);
                        % st2 = OTFS_modulation(N,M,(X));
                    end
                end
                X_td = X.*(twiddle_factor);
                X_tw = twisted_conv(X_td, g_tx, N, M);    % 发射端：扭曲卷积 + IDZT
                % s = OTFS_modulation(N,M,(X));
                s_mat=ifft(X_tw, [],2)*sqrt(N);
                s = s_mat(:);

                if num3==2
                    if  cc>=2
                        st2=s;
                    end
                end

                %% OTFS channel output%%%%%
                r = OTFS_channel_output(N,M,taps,delay_taps,Doppler_taps,chan_coef,s,SNR_dB(num1),cpl);
                for num4=1:taps
                    Y(num4,(num3-1)*(M*N+cpl)+1+delay_taps(num4):(num3-1)*(M*N+cpl)+N*M+delay_taps(num4)+cpl)=r(num4,1+delay_taps(num4):N*M+delay_taps(num4)+cpl);
                end
            end
            L_s=randi(20);% Remove part of the first frame

            % Z=Y(:,1:3*(M*N+cpl)+cpl);
            Z=Y;
            R=sum(Z,1);
            R_temp=awgn(R,SNR_dB(num1));   %Noise
            % R_temp=R;
            R=R_temp(L_s+1:end);

            sx=st1;
            if  cc>=2
                sx=st2;
            end
            
            Lt=max(delay_taps);
            Kv=ceil(max(k_max));
            % doppler_taps_set=[-Kv:0.1:Kv].';     % 正常来说应该是这个，但是生成的信道的多普勒没有负数
            rv=0.05;                % 一般来说值越小，同步性能越好，但复杂度越高
            doppler_taps_set=[-Kv:rv:Kv].';

            A_pccf=zeros(length(R),length(doppler_taps_set));
            sx1 = [sx(N*M-cpl+1:N*M);sx];%add one cp
            for ii=1:length(doppler_taps_set)

                doppler_taps_i=doppler_taps_set(ii);
                % sx_=(sx.').*exp(1j*2*pi/M*(0:length(sx)-1)*doppler_taps_i/N);
                sx_=(sx1.').*exp(1j*2*pi/M*(-cpl:-cpl+length(sx1)-1)*doppler_taps_i/N);
                sx_temp=[sx_,zeros(1,length(R)-length(sx_))];
                KK2_=pccf(sx_temp,R); %Matched Filter , calculating PCCF
                A_pccf(:,ii)=KK2_;
            end
    
            ind_true=M*N+2*cpl-L_s+1;

           

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [max_1,idx]=max(abs(A_pccf));
            [max_2,idx_doppler]=max(abs(max_1));
          %  [maxValue, linearIndex] = max(data(:));  % 先找到线性索引
            row_ind = idx(idx_doppler);
            KK2=A_pccf(:,idx_doppler);% 这里是补偿后的，第二次尝试可以是未补长的 ，已经得到重点row_ind 位置
            % plot(abs(KK2))
            

            % F=abs(KK2(1:end-M*N));
            F=abs(KK2);
            %Find possible locations
            [sortedValues, sortedIndices] = sort(abs(F), 'descend');
            ind_max=sortedIndices(1);
            if ind_max<=7
                ind_min=1;
            else
                ind_min=ind_max-Lt;
            end

            if(KP(cc,2)==1)%Find the first 'bignum' largest numbers
                bignum=KP(cc,1);
            else
                bignum=KP(cc,1);
            end

            if(KP(cc,1)==4)
                bignum=5;
            end


            top5Values = sortedValues(1:bignum);
            top5Indices = sortedIndices(1:bignum);


            Ksum=zeros(length(top5Indices),1);
            for numfind=1:length(top5Values)
                if top5Indices(numfind)<11
                    continue;
                end


                Ksum(numfind)=sum(F(top5Indices(numfind)-8:top5Indices(numfind)-1));
            end

            A = Ksum;
            nonZeroValues = A(A ~= 0);

            for kk=1:length(nonZeroValues)
                [minNonZeroValue,pos_min] = min(nonZeroValues);
                minNonZeroIndex = find(A == minNonZeroValue);

                ind_temp=top5Indices(minNonZeroIndex);
                if ind_temp >= ind_min && ind_temp <= ind_max
                    break
                else
                    nonZeroValues(pos_min)=NaN;
                end
            end


            %Determine location
            %% 完美同步的BER

            r_1=sum(r);
            % noise_1 = sqrt(sigma2/2)*(randn(size(r_1)) + 1i*randn(size(r_1)));
            % r_1=r_1+noise_1;
            r_1=awgn(r_1,SNR_dB(num1));   %Noise
            r_2=r_1(cpl+1:cpl+(N*M));
            if cc==1
                r_mat = reshape(r_2,M,N);  % 串并转换
                Y_DD = fft(r_mat, [],2)/sqrt(N);
                R_matched = twisted_conv(Y_DD, conj(g_rx), N, M) .* conj(twiddle_factor);
                y=R_matched(:);

                
                x_est=Wmmse_tradition*y;
                % x_est = OTFS_mp_detector(N,M,M_mod,taps,delay_taps,Doppler_taps,chan_coef,sigma2,Y);
                % [~,x_est] = MPA_detector(N,M,M_mod,sigma2,y,H,50);
                data_demapping = qamdemod(x_est,M_mod,'gray');
                data_info_est = reshape(de2bi(data_demapping,M_bits),N_bits_perfram,1);
                errors = sum(xor(data_info_est,data_info_bit));
                ber_true=errors/N_bits_perfram;
                ber_fram_true(num2)=ber_true;
            end

            %% 非完美同步的BER
            %%%%%%%%
            % ind_true
            ind=top5Indices(minNonZeroIndex);
            ind=ind+cpl;


            if ind==ind_true
                %%% 正确同步
                KK=KK+1;
            else
                %%% 未正确同步
                is_syn(num2,1)=1;
                is_syn1(num2,cc)=1;
                if ind+N*M+cpl+N*M-1>length(R)    % 为了保证能得到一个完整的OTFS帧
                    ind =min(top5Indices);
                    if ind+N*M+cpl+N*M-1>length(R)
                        ind=floor(ind/2);
                    end
                end

                pos_syn=[pos_syn; ind-ind_true];
            end

            %%%%%%%%%%%%%%%%%%%%%% 计算BER
            r_temp=R(ind+N*M+cpl:ind+N*M+cpl+N*M-1);
            % r_temp=r_2;
            r_mat = reshape(r_temp,M,N);  % 串并转换
            Y_DD = fft(r_mat, [],2)/sqrt(N);
            R_matched = twisted_conv(Y_DD, conj(g_rx), N, M);
            R_matched=(R_matched).*conj(twiddle_factor);
            y=R_matched(:);
            % x_est = OTFS_mp_detector(N,M,M_mod,taps,delay_taps,Doppler_taps,chan_coef,sigma2,Y);
            % [~,x_est] = MPA_detector(N,M,M_mod,sigma2,y,H,50);
            % Wmmse_tradition=(HDD'*HDD+sigma2*eye(N*M))\HDD';
            x_est=Wmmse_tradition*y;
            data_demapping = qamdemod(x_est,M_mod,'gray');
            data_info_est = reshape(de2bi(data_demapping,M_bits),N_bits_perfram,1);
            errors = sum(xor(data_info_est,data_info_bit));
            ber=errors/N_bits_perfram;
            % if cc==1
            %     ber_snr(num2)=ber_true;
            %     ber_snr_1(num2,cc)=ber_true;
            % else
            %     ber_snr(num2)=ber;
            %     ber_snr_1(num2,cc)=ber;
            % end
            ber_snr(num2,1)=ber;
            ber_snr_1(num2,cc)=ber;


            
            fprintf('OTFS: %d   %d   %d\n',cc,num1,num2);


        end

        To(cc,num1)=KK/BL;

        if cc==1
            [idx,~]=find(is_syn(:,1));
            ber_snr_no_syn_err=ber_snr;
            ber_snr_no_syn_err(idx)=[];
            BER_no_syn_err(num1)=mean(ber_snr_no_syn_err);
        end

        BER(num1,cc)=mean(ber_snr);

        if cc==1
            BER_true(num1,cc)=mean(ber_fram_true);
        end

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
semilogy(SNR_dB,BER_true,'-*g');
hold on;
semilogy(SNR_dB,BER(:,2),'-*b');
semilogy(SNR_dB,BER(:,1),'-*r');
grid on;
xlabel('SNR (dB)');
ylabel('BER');
legend("Perfect synchronization",'Random Zak sequence','Proposed Zak sequence');
filename1=strcat('ZAK-OTFS仿真结果\',string(datetime('now','Format','yyyy-MM-dd''T''HHmmss')),'_BER-cp-3-1.fig');
saveas(gcf,filename1)

figure(2)
semilogy(SNR_dB,To(1,:),'-*r');
hold on;
semilogy(SNR_dB,To(2,:),'-*b')
ax = gca;
ax.XTick = SNR_dB(1):5:SNR_dB(end);
ax.YTick = 0:0.2:1;
ax.FontName = 'Times New Roman'; % 字体类型
ax.YTickLabel = {'0%','20%','40%','60%','80%','100%'};
ax.YScale="linear";
ax.XScale="linear";
grid on;
ylim([0 1]);
xlim([SNR_dB(1) SNR_dB(end)]);
xlabel('SNR/dB');
ylabel('Sync success prob');
legend('Proposed Zak sequence', 'Random Zak sequence');
filename1=strcat('ZAK-OTFS仿真结果\',string(datetime('now','Format','yyyy-MM-dd''T''HHmmss')),'_syn-cp-3-1.fig');
saveas(gcf,filename1)





%% -----------------------
% 扭曲卷积函数
function Y = twisted_conv(X, g, N, M)
Y = zeros(M,N);
for n = 0:N-1
    for m = 0:M-1
        tmp = 0;
        for np = 0:N-1
            for mp = 0:M-1
                nn = mod(n-np,N);
                mm = mod(m-mp,M);
                tmp = tmp + X(mp+1,np+1) * g(mm+1,nn+1) * exp(1j*2*pi*nn*mp/M);
            end
        end
        Y(m+1,n+1) = tmp;
    end
end
end


function r = OTFS_channel_output0(N,M,taps,delay_taps,Doppler_taps,chan_coef,sigma_2,s)
%% wireless channel and noise
L = max(delay_taps);
s = [s(N*M-L+1:N*M);s];%add one cp
s_chan = 0;
for itao = 1:taps
    s_chan = s_chan+chan_coef(itao)*circshift([s.*exp(1j*2*pi/M ...
        *(-L:-L+length(s)-1)*Doppler_taps(itao)/N).';zeros(delay_taps(end),1)],delay_taps(itao));
end
noise = sqrt(sigma_2/2)*(randn(size(s_chan)) + 1i*randn(size(s_chan)));
r = s_chan + noise;
r = r(L+1:L+(N*M));%discard cp
end


function H_eff = compute_H_eff_ZAK_OTFS(N, M, delays,dopplers,gains,P, g_tx, g_rx,twiddle_factor,numSym,I_F_inv,I_F)
% H_eff = zeros(numSym,numSym);
% for k = 1:numSym
%     e = zeros(numSym,1); e(k)=1;
%     % 发射端卷积 + IDZT
%     E=reshape(e,N,M).*(twiddle_factor.');
%     Y_E=ifft(twisted_conv(E, g_tx,N,M), N,1)*sqrt(N);
%     s_e = reshape(Y_E.',[],1);
%     r_e = OTFS_channel_output(N,M,P,delays,dopplers,gains,0,s_e);
%     % 接收端匹配卷积
%     r_e_blocks = reshape(r_e,M,N);
%     R_e = fft(r_e_blocks.',N,1)/sqrt(N);
%     R_e_matched = twisted_conv(R_e, conj(g_rx), N, M) .* conj(twiddle_factor.');
%     H_eff(:,k) = R_e_matched(:);
% end

tw_vec = twiddle_factor.'; % N x M
tw_vec = tw_vec(:);        % NM x 1
D_tw_tx = spdiags(tw_vec,0,numSym,numSym);       % 发送端乘 twiddle
D_tw_rx = spdiags(conj(tw_vec),0,numSym,numSym); % 接收端乘 conj(twiddle)

%  构造扭曲卷积矩阵 T_tx 和 T_rx，使 vec(Y)=T*vec(X)
T_tx = build_twisted_conv_mat(g_tx, N, M);          % NM x NM
T_rx = build_twisted_conv_mat(conj(g_rx), N, M);    % NM x NM

% 置换矩阵 K such that vec(A') = K * vec(A)
perm_idx = reshape(1:numSym, N, M).'; % M x N then transpose-> to get mapping
perm_idx = perm_idx(:);
K = sparse(1:numSym, perm_idx, 1, numSym, numSym); % vec(A') = K * vec(A)

% B: 从 vec(X) (DD-domain) 到时域发送 s 的线性映射
% steps: vec(X) -> multiply twiddle -> twisted_conv (T_tx) -> apply IFFT on rows -> reshape to time vector s
B = K * I_F_inv * T_tx * D_tw_tx; % size: (NM x NM) mapping vec(X) -> s (length NM)

% A: 从时域接收 r (length NM) 到 vec(R_matched) 的线性映射
% steps: r -> reshape -> transpose -> FFT on rows -> twisted_conv (T_rx) -> multiply conj(twiddle)
A = D_tw_rx * T_rx * I_F * K; % size: (NM x NM), mapping r -> vec(R_matched)

%  构造时域信道矩阵 C (txLen x txLen)
% 通过对每个时域单位脉冲输入调用 OTFS_channel_output 构造 C 的每一列
H_T=zeros(M*N, M*N);
for i=1:P
    T=OTFS_Tmartix_gen(N,M,delays,dopplers,i);
    H_T=H_T+gains(i)*T;
end
H_eff = A * H_T * B;  % NM x NM
end



function T=OTFS_Tmartix_gen(N,M,delay_taps,Doppler_taps,i)

%生成第i径的T矩阵
%N M为延迟多普勒域维度
%delay_taps为延迟时间
%Doppler_taps为多普勒频移
%i 为第几径
%% 首先生成单位阵 I_m I_n
% Im=eye(M);
% In=eye(N);
%
% %% 生成DFT矩阵
% Fn=exp(-1i*2*pi./N*((0:N-1).'*(0:N-1)));

%% 生成Delay矩阵 前向循环移位
Temp_permutation_matrix=eye(N*M);
Temp_permutation_matrix=circshift(Temp_permutation_matrix(1:N*M,:),delay_taps(i));

%% 生成Doppler矩阵  对角阵时变系数
delta_matrix=diag(exp(1i*2*pi*Doppler_taps(i)*(0:M*N-1)./N./M));

%% 按照公式依次生成P Q T矩阵
%% 生成P矩阵
P=Temp_permutation_matrix;

%% 生成Q矩阵
Q=delta_matrix;

T=P*Q;

end


function T = build_twisted_conv_mat(g, N, M)
% g: N x M (或按使用的索引)
NM = N*M;
T = zeros(NM, NM);
% 索引约定：把 N x M 的矩阵按 MATLAB 的 column-major 展开为 vec(X),
% 元素 (n+1, m+1) 的线性索引为 n+1 + m*N (n from 0..N-1, m from 0..M-1)
for n = 0:N-1
    for m = 0:M-1
        out_idx = n + 1 + m*N;
        for np = 0:N-1
            for mp = 0:M-1
                in_idx = np + 1 + mp*N;
                nn = mod(n - np, N);
                mm = mod(m - mp, M);
                coeff = g(nn+1, mm+1) * exp(1j*2*pi*nn*mp/M);
                T(out_idx, in_idx) = coeff;
            end
        end
    end
end
end
