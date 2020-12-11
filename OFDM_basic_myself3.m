%%%%%%%%%%%%%%%%%%%%%       ��������źż���ཻ��    %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%        OFDM_basic_myself3.m            %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%      data:2020��12��11��  author:����󽫾� %%%%%%%%%%

%%%%%%%%%%%%%%%%%����˵��
%%%��1��CP��ZP����   ��2�����ز���������ż���ŵ���������������

%%%%%%    ���滷��
%�����汾��MATLAB R2019a

%**************************** �������� **************************%
clear all;

%%%%%%%%%%%%%%%%%%%%�����趨%%%%%%%%%%%%%


%%%%ѡ��CP��ZP
NgType = 1;  %����ZP��CP NgType = 1��2
if NgType == 1
    nt = 'CP';
elseif NgType == 2
    nt = 'ZP';
end

%%%%ѡ���ŵ�����
Ch = 0;
if Ch == 0
    chType ='AWGN'; %��˹�������ŵ� 
    Target_neb = 100;
else
    chType ='CH';
    Target_neb = 500;
end
figure(Ch+1);
clf;

PowerdB = [0 -8 -17 -21 -25]; %�ŵ���ͷ��������'dB'
Delay = [0 3 5 6 8]; %�ŵ�ʱ�� 
Power = 10.^(PowerdB/10); %�ŵ���ͷ�������� '����'
Ntap = length(PowerdB);
Lch = Delay(end)+1;
Nbps = 2;   %���ƽ��� 2/4/6
M = 2^Nbps; %QPSK��16-QAM��64-QAM
Nfft = 64; %FFT��С

%Ng = 3;

Ng = Nfft/4; %���������GI�����ȣ���û�б��������Ng = 0
Nsym = Nfft + Ng; %��������

%%%����Nvc
Nvc = Nfft/4; %Nvc������0����û��VC���������ز���
% Nvc = 0; %Nvc������0����û��VC���������ز���

Nused = Nfft - Nvc; %NusedΪ���ڴ������ݵ����ز���

EbN0 = [0:2:20];  %Eb/N0
% EbN0 = 50;  %Eb/N0
N_iter = 1e5;  %����ÿһ��EbN0�ĵ�������
Nframe = 3; %ÿһ֡�ķ�����
sigPow = 0; %��ʼ�źŹ���
file_name = ['OFDM_BER_' chType '_' nt '_' 'GL' num2str(Ng) '.dat'];
fid = fopen(file_name,'w+');
norms = [1 sqrt(2) 0 sqrt(10) 0 sqrt(42)];  %BPSK 4-QAM 16-QAM


for i = 0:length(EbN0)
    randn('state',0);
    rand('state',0);
    Ber2=ber(); %��ʼ��BER
    Neb = 0; %��ʼ���������
    Ntb = 0; %��ʼ���ܱ�����
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
        
        %%%%%�������
        if i == 0   %ֻ�����źŹ���
            sigPow_temp = x_GI*x_GI';;
        end
        
        
        
        if Ch==0
            y = x_GI;  %û���ŵ�
        else      %�ྶ˥���ŵ�
            channel =(randn(1,Ntap)+1j*randn(1,Ntap)).*sqrt(Power/2);
            h = zeros(1,Lch);
            h(Delay+1) = channel;
            y = conv(x_GI,h);
        end
        
        if i == 0   %ֻ�����źŹ���
            y1 = y(1:Nframe*Nsym);
            sigPow = sigPow + y1*y1';
            continue;
        end

        %******************** �ŵ� ***********************%
        snr = EbN0(i) + 10*log10(Nbps*(Nused/Nfft));  %%����fig��ţ���1��,ԭ�鹫ʽ
%         snr = EbN0(i) + 10*log10(Nbps);
%         snr = EbN0(i) + 10*log10(Nbps*(Nfft/Nsym));  %%����fig��ţ���3��CP��������
        noise_msg = sqrt((10.^(-snr/10))*sigPow/2);
        y_GI = y + noise_msg*(randn(size(y)) + 1j*randn(size(y)));
        
        %%%%%%%%%%%���ն�
        kk1 = (NgType==2)*Ng + [1:Nsym];
        kk2 = 1:Nfft;
        kk3 = 1:Nused;
        kk4 = Nused/2 + Nvc + 1:Nfft;
        kk5 = (Nvc~=0)+[1:Nused/2];
        if Ch ==1
            H = fft([h zeros(1,Nfft-Lch)]);  %�ŵ�Ƶ����Ӧ
            H_shift(kk3) = [H(kk4) H(kk5)];
        end
        
        for k =1:Nframe
            Y(kk2) = fft(remove_GI(Ng,Nsym,NgType,y_GI(kk1)));
            Y_shift = [Y(kk4) Y(kk5)];
            if Ch ==0
                Xmod_r(kk3) = Y_shift;
            else
                Xmod_r(kk3) = Y_shift./H_shift;  %������
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

%%%%ѡ��CP��ZP
if Ch == 1
    if NgType == 1
        a = load(file_name);
        save('ofdm_basic_myself3_cp16_rayleigh','a');
    elseif NgType == 2
        b = load(file_name);
        save('ofdm_basic_myself3_zp16_rayleigh','b');
    end
else
    

%%%%%%%%%%%%%%ʵ���¼
%%%%2020��12��11��
%%%%CP��ZP�����𵽱�����������ã����������ز����źͷ��ż����
%%%%��AWGN�������ŵ��¾��ܻ�����ȷ����������

    

        
            
   