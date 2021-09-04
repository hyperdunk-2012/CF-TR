%% QPSK
clc
clear
close all

num = 100000;
K = 64;
N = 58;
L = 4;

filename = 'test.mat';
QPSK_SET = [1+1i,1-1i,-1+1i,-1-1i];

pos1 = randperm(K,K);
pos0 = zeros(1,K-2);
m = 1;
for n = 1:K
    if pos1(n) ~= 1 && pos1(n) ~= (K/2)+1
        pos0(m) = pos1(n);
        m=m+1;
    end
end
pos = pos0(1:N);
TR_pos = pos0(N+1:N+L);
pos = sort(pos);
TR_pos = sort(TR_pos);

SymbolTX = zeros(K,num);

for n = 1:num
    for m = 1:N
        SymbolTX(pos(m),n) = QPSK_SET(randi([1,4]));
    end    
end

TX = ifft(SymbolTX,2*K,1);

Signal_Power = abs(TX.^2);
Peak_Power   = max(Signal_Power,[],1);
Mean_Power   = mean(Signal_Power,1);
PAPR_Original = 10*log10(Peak_Power./Mean_Power);

[cdf1, PAPR1] = ecdf(PAPR_Original);
ecdf1 = 1-cdf1;

figure(1)
semilogy(PAPR1,ecdf1,'-b')
ylim([0.001,1]);
legend('Original')
title('ORIGIN')
xlabel('PAPR0 [dB]');
ylabel('CCDF (Pr[PAPR>PAPR0])');
grid on

save(filename);