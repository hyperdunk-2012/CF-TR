%% CF-TR
I = 3;
% M = 5.3;

SymbolTR = zeros(K,num);
SymbolTRAIN = zeros(K,num);

tic

for m = 1:num
    S = SymbolTX(1:K,m);
    s = TX(1:K,m);
    f = zeros(K,1);
    T = zeros(K,1);
    s0 = zeros(K,1);
    Signal_Power = abs(s.^2); 
    Mean_Power   = mean(Signal_Power);
    
    for n = 1:I

        M = sqrt(Mean_Power*log(K/L));
        
        for k = 1:K
            if abs(s(k,1)) >= M
                s0(k,1) = M*exp(1i*phase(s(k,1)));
            else
                s0(k,1) = s(k,1);
            end
        end
        
        f = s0-s;
        
        k0 = 1;
        clear f0;
        clear pos0;
        
        for k = 1:K
            if f(k,1) ~= 0
                pos0(1,k0) = k;
                f0(k0,1) = f(k,1);
                k0 = k0+1;
            end
        end
        
        k0 = k0-1;
        
        if k0 ~= 0
            R = zeros(k0,L);

            for k = 1:k0
                for k1 = 1:L 
                    R(k,k1) = exp(1i*2*pi*(TR_pos(k1)-1)*(pos0(k)-1)/K)/sqrt(K);
                end
            end

            R0 = pinv(R);
            F = R0*f0;

            for k = 1:L
                S(TR_pos(k),1) = F(k,1);
                T(TR_pos(k),1) = F(k,1);
            end
        end
        
        s = sqrt(K)*ifft(S);

    end
    
    SymbolTR(1:K,m) = S;
    SymbolTRAIN(1:K,m) = T;
    
end

toc
tim=toc;

TRAIN = sqrt(K)*ifft(SymbolTRAIN,K,1);
TR = sqrt(K)*ifft(SymbolTR,K,1);

Signal_Power = abs(TR.^2);
Peak_Power   = max(Signal_Power,[],1);
Mean_Power   = mean(Signal_Power,1);
PAPR_TR2 = 10*log10(Peak_Power./Mean_Power);

[cdf2, PAPR2] = ecdf(PAPR_TR2);
ecdf2 = 1-cdf2;

figure(2)
semilogy(PAPR1,ecdf1,'-b',PAPR2,ecdf2,'-r')
ylim([0.001,1]);
legend('Original','CF-TR')
title('COMPARE')
xlabel('PAPR0 [dB]');
ylabel('CCDF (Pr[PAPR>PAPR0])');
grid on

save(filename);