clc
clearvars
close all

load train;
unquantizedSignal = y;
mean = mean(y);
Vp = 1;
N = 4;


[quantizedSignal,varLin,varSat,varTeo,SNqR,SNqRTeo] = MyQuantizer(unquantizedSignal,Vp,N);
estimatedBitStream = MyGraycode(quantizedSignal,Vp,N);
estimatedSignal = MyDAconverter(estimatedBitStream,Vp,N);

sound(estimatedSignal)

%results
figure
plot(unquantizedSignal)
title("Unquantized Signal")
ylabel("Amplitude (V)")
xlabel("Time step (1/Fs)")

figure
plot(quantizedSignal)
title(['Quantized Signal, Vp = ',num2str(Vp),', N = ',num2str(N)])
ylabel("Amplitude (V)")
xlabel("Time step (1/Fs)")

figure
plot(estimatedSignal)
title(['Estimated Signal, Vp = ',num2str(Vp),', N = ',num2str(N)])
ylabel("Amplitude (V)")
xlabel("Time step (1/Fs)")

x = 6000:6041;
figure
plot(x,unquantizedSignal(6000:6041))
hold on
stairs(x,estimatedSignal(6000:6041))
xlim([6000,6040])
ylim([-Vp,Vp])
title(['Unquantized vs Estimated Signal, Vp = ',num2str(Vp),', N = ',num2str(N)])
legend("Unquantized Signal","Estimated Signal")
legend('Location','northwest')
ylabel("Amplitude (V)")
xlabel("Time step (1/Fs)")

% figure
% plot(SNqR)
% title(['Signal to Noise Ratio, Vp = ',num2str(Vp),', N = ',num2str(N)])
% ylabel("Amplitude (dB)")
% xlabel("Time step (1/Fs)")

function [quantizedSignal,varLin,varSat,varTeo,SNqR,SNqRTeo] = MyQuantizer(unquantizedSignal,Vp,N)
    %quantization setup
    levels = 2^N;
    step = 2*Vp/levels;
    varTeo = step^2/12;
    SNqRTeo = mag2db(levels^2);
    level = linspace(-(levels/2-0.5)*step,(levels/2-0.5)*step,levels);
    quantizedSignal = nan(1,length(unquantizedSignal));
    
    %quantization
    for i = 1:length(unquantizedSignal)
        for j = 1:length(level)
            if unquantizedSignal(i) <= level(j)+step/2 && unquantizedSignal(i) > level(j)-step/2
                quantizedSignal(i) = level(j);
            end
        end  
        if unquantizedSignal(i) < level(1)
                quantizedSignal(i) = level(1);
        elseif unquantizedSignal(i) > level(length(level))
                quantizedSignal(i) = level(length(level));     
        end
    end
    
    %Variance linear
    varLin=var(quantizedSignal.' - unquantizedSignal);

    % Saturated error variance
    satError = quantizedSignal.' - unquantizedSignal;
    satError = min(max(satError, -Vp), Vp);
    varSat = var(satError);

    % Signal to Quantization Noise power Ratio (SNqR) in dB
    SNqR = 20 * log10(var(unquantizedSignal) ./ varSat);
end

function [bitStream] = MyGraycode(quantizedSignal,Vp,N)
    levels = 2^N;
    step = 2*Vp/levels;
    bitStream = nan(1,N*length(quantizedSignal));
    for i = 1:length(quantizedSignal)
        bit = (quantizedSignal(i)+(levels/2-0.5)*step)/step;
        bin = dec2bin(bit,N);
        bitStream(1,1+N*(i-1)) = str2double(bin(1));
        for j = 2:N
            bitStream(1,j+N*(i-1)) = bitxor(str2double(bin(1,j)),str2double(bin(1,j-1)));
        end    
    end
end

function [estimatedSignal] = MyDAconverter(estimatedBitStream,Vp,N)
    levels = 2^N;
    step = 2*Vp/levels;
    level = linspace(-(levels/2-0.5)*step,(levels/2-0.5)*step,levels);
    signalbits = nan(1,length(estimatedBitStream));
    estimatedSignal = zeros(1,length(estimatedBitStream)/N);
    for i = 1:length(estimatedSignal)
        signalbits(1+N*(i-1)) = estimatedBitStream(1+N*(i-1));
        for j = 2:N
            signalbits(j+N*(i-1)) = bitxor(estimatedBitStream(j+N*(i-1)),signalbits(j+N*(i-1)-1));
        end
        for k = 1:N
             estimatedSignal(i) = estimatedSignal(i) + signalbits(k+N*(i-1))*2^(N-k);
        end
    end
    for l = 1:length(estimatedSignal)
        estimatedSignal(l) = level(estimatedSignal(l)+1);
    end 
end



