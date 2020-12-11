%%%%%%%%%%%%%%%%%%%%%       检查正弦信号间的相交性    %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%        OFDM_basic_myself3.m            %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%      data:2020年12月11日  author:飞蓬大将军 %%%%%%%%%%

%%%%%%%%%%%%%%%%%程序说明
%%%（1）CP和ZP作用   （2）子载波干扰与符号间干扰的区别及其消除方法

%%%%%%    仿真环境
%软件版本：MATLAB R2019a

%**************************** 程序主体 **************************%
clear all;

%%%%%%%%%%%%%%%%%%%%参数设定%%%%%%%%%%%%%


%%%%选择CP或ZP
NgType = 1;  %对于ZP或CP NgType = 1或2
if NgType == 1
    nt = 'CP';
elseif NgType == 2
    nt = 'ZP';
end

%%%%选择信道类型
Ch = 0;
if Ch == 0
    chType ='AWGN'; %高斯白噪声信道 
    Target_neb = 100;
else
    chType ='CH';
    Target_neb = 500;
end
figure(Ch+1);
clf;

PowerdB = [0 -8 -17 -21 -25]; %信道抽头功率特性'dB'
Delay = [0 3 5 6 8]; %信道时延 
Power = 10.^(PowerdB/10); %信道抽头功率特性 '线性'
Ntap = length(PowerdB);
Lch = Delay(end)+1;
Nbps = 2;   %调制阶数 2/4/6
M = 2^Nbps; %QPSK、16-QAM、64-QAM
Nfft = 64; %FFT大小

%Ng = 3;

Ng = Nfft/4; %保护间隔（GI）长度，若没有保护间隔，Ng = 0
Nsym = Nfft + Ng; %符号周期

%%%调整Nvc
Nvc = Nfft/4; %Nvc若等于0，则没有VC（虚拟子载波）
% Nvc = 0; %Nvc若等于0，则没有VC（虚拟子载波）

Nused = Nfft - Nvc; %Nused为用于传输数据的子载波数

EbN0 = [0:2:20];  %Eb/N0
% EbN0 = 50;  %Eb/N0
N_iter = 1e5;  %对于每一次EbN0的迭代次数
Nframe = 3; %每一帧的符号数
sigPow = 0; %初始信号功率
file_name = ['OFDM_BER_' chType '_' nt '_' 'GL' num2str(Ng) '.dat'];
fid = fopen(file_name,'w+');
norms = [1 sqrt(2) 0 sqrt(10) 0 sqrt(42)];  %BPSK 4-QAM 16-QAM


for i = 0:length(EbN0)
    randn('state',0);
    rand('state',0);
    Ber2=ber(); %初始化BER
    Neb = 0; %初始化误比特数
    Ntb = 0; %初始化总比特数
    for m = 1:N_iter
        X = randi([0 M-1],1,Nused*Nframe);
        Xmod = qammod(X,M,'gray')/norms(Nbps);
        if NgType~=2
            x_GI = zeros(1,Nframe*Nsym);
        elseif  NgType == 2
            x_GI = zeros(1,Nframe*Nsym+Ng);
        end
%         kk1 = 1:Nused/2;
%         kk2 = Nused/2+1:Nused;
        kk1 = [1:Nused/2];
        kk2 = [Nused/2+1:Nused];
        kk3 = 1:Nfft;
        kk4 = 1:Nsym;
        for k = 1:Nframe
            if Nvc~= 0
                X_shift = [0 Xmod(kk2) zeros(1,Nvc-1) Xmod(kk1)];
            else
                X_shift = [Xmod(kk2) Xmod(kk1)];
            end
            x = ifft(X_shift);
            x_GI(kk4) = guard_interval(Ng,Nfft,NgType,x);
            kk1 = kk1 + Nused;
            kk2 = kk2 + Nused;
            kk3 = kk3 + Nfft;
            kk4 = kk4 + Nsym;
        end
        
        %%%%%能量检测
        if i == 0   %只测量信号功率
            sigPow_temp = x_GI*x_GI';;
        end
        
        
        
        if Ch==0
            y = x_GI;  %没有信道
        else      %多径衰落信道
            channel =(randn(1,Ntap)+1j*randn(1,Ntap)).*sqrt(Power/2);
            h = zeros(1,Lch);
            h(Delay+1) = channel;
            y = conv(x_GI,h);
        end
        
        if i == 0   %只测量信号功率
            y1 = y(1:Nframe*Nsym);
            sigPow = sigPow + y1*y1';
            continue;
        end

        %******************** 信道 ***********************%
        snr = EbN0(i) + 10*log10(Nbps*(Nused/Nfft));  %%方便fig标号，（1）,原书公式
%         snr = EbN0(i) + 10*log10(Nbps);
%         snr = EbN0(i) + 10*log10(Nbps*(Nfft/Nsym));  %%方便fig标号，（3）CP消耗能量
        noise_msg = sqrt((10.^(-snr/10))*sigPow/2);
        y_GI = y + noise_msg*(randn(size(y)) + 1j*randn(size(y)));
        
        %%%%%%%%%%%接收端
        kk1 = (NgType==2)*Ng + [1:Nsym];
        kk2 = 1:Nfft;
        kk3 = 1:Nused;
        kk4 = Nused/2 + Nvc + 1:Nfft;
        kk5 = (Nvc~=0)+[1:Nused/2];
        if Ch ==1
            H = fft([h zeros(1,Nfft-Lch)]);  %信道频率响应
            H_shift(kk3) = [H(kk4) H(kk5)];
        end
        
        for k =1:Nframe
            Y(kk2) = fft(remove_GI(Ng,Nsym,NgType,y_GI(kk1)));
            Y_shift = [Y(kk4) Y(kk5)];
            if Ch ==0
                Xmod_r(kk3) = Y_shift;
            else
                Xmod_r(kk3) = Y_shift./H_shift;  %均衡器
            end
            kk1 = kk1 + Nsym;
            kk2 = kk2 + Nfft;
            kk3 = kk3 + Nused;
            kk4 = kk4 + Nfft;
            kk5 = kk5 + Nfft;
        end
        X_r = qamdemod(Xmod_r*norms(Nbps),M,'gray');
        Neb = Neb + sum(sum(de2bi(X_r,Nbps)~=de2bi(X,Nbps)));
        Ntb = Ntb + Nused*Nframe*Nbps;
%         if Neb>Target_neb
%             break
%         end
    end
    if i == 0
        sigPow = sigPow/Nsym/Nframe/N_iter;
        fprintf('Signal power= %11.3e\n', sigPow);
        fprintf(fid,'%%Signal power= %11.3e\n%%EbN0[dB]       BER\n', sigPow);
    else
        Ber = Neb/Ntb;
        fprintf('EbN0=%3d[dB], BER=%4d/%8d=%11.3e\n', EbN0(i), Neb,Ntb,Ber)
        fprintf(fid, '%d\t%11.3e\n', EbN0(i), Ber);
        if Ber<1e-6
            break;
        end
    end
end %end for i
if(fid~=0)
     fclose(fid);
end
disp('sumualtion is finished');
plot_ber(file_name,Nbps);

%%%%选择CP或ZP
if Ch == 1
    if NgType == 1
        a = load(file_name);
        save('ofdm_basic_myself3_cp16_rayleigh','a');
    elseif NgType == 2
        b = load(file_name);
        save('ofdm_basic_myself3_zp16_rayleigh','b');
    end
else
    

%%%%%%%%%%%%%%实验记录
%%%%2020年12月11日
%%%%CP和ZP都能起到保护间隔的作用，消除了子载波干扰和符号间干扰
%%%%在AWGN和瑞利信道下均能画出正确误码率曲线

    

        
            
   