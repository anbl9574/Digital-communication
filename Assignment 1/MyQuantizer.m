function[quantizedSignal,varLin,varSat,SNqR] = MyQuantizer(unquantizedSignal,Vp,N)

step=(Vp*2)/2^N;
maxlvl = (2^N/2-0.5)*step;
minlvl = -(2^N/2-0.5)*step;

levelquant = linspace(minlvl,maxlvl,2^N);

holder=zeros(length(unquantizedSignal));

for i=1:length(unquantizedSignal)

    for j=1:length(levelquant)

        holder(i,j)=abs(unquantizedSignal(i)-levelquant(j));

    end 
end 


[~, minIndices] = min(holder, [], 2);

quantizedSignal=levelquant(minIndices);

%Variance linear

varLin=var(quantizedSignal-unquantizedSignal);

% Saturated error variance
satError = quantizedSignal - unquantizedSignal;
satError = min(max(satError, -Vp), Vp);
varSat = var(satError);

% Signal to Quantization Noise power Ratio (SNqR) in dB
SNqR = 10 * log10(var(unquantizedSignal) / varSat);

t = linspace(0, 1, numel(unquantizedSignal));
figure;
subplot(2, 1, 1);
plot(t, unquantizedSignal, '-b', t, quantizedSignal, '--r');
title('Original and Quantized Signals');
legend('Original Signal', 'Quantized Signal');
xlabel('Time');
ylabel('Amplitude');