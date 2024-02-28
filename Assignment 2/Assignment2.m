clc
clearvars
close all

load train;
unquantizedSignal = y;
matchedFilterFlag = 1;
Vp = 1;
N = 3;
Es = 1;
M = 2^N;

% From Assignment 1
[quantizedSignal,~,~,~,~,~] = MyQuantizer(unquantizedSignal,Vp,N);
transmitedBitStream = MyGraycode(quantizedSignal,Vp,N);

% Assignment 2
transmitSignal = MyMPAM(transmitedBitStream,M,Es);
receiveSignal = MyAWGNchannel(transmitSignal,0.01);
[estimatedBitStream,BER1] = DemodulateMPAM(receiveSignal, M, Es, transmitedBitStream, matchedFilterFlag);
estimatedSignal = MyDAconverter(estimatedBitStream,Vp,N);

% Random index for plots
sx = randi(length(y) -1000);

figure
% Plot train signal
plot(unquantizedSignal)
title("Unquantized Signal")
ylabel("Amplitude (V)")
xlabel("Time step (1/Fs)")

figure
% Plot quantized train signal
plot(quantizedSignal)
title(['Quantized Signal, Vp = ',num2str(Vp),', N = ',num2str(N)])
ylabel("Amplitude (V)")
xlabel("Time step (1/Fs)")

figure;
% Plot received signal with AWGN
plot(receiveSignal(sx:sx+10000), 'LineWidth', 2, 'DisplayName', 'Received Signal with AWGN');
hold on;
% Plot transmitted signal without AWGN with transparency
plot(transmitSignal(sx:sx+10000), 'LineWidth', 2, 'DisplayName', 'Transmitted Signal', 'Color',[1 0.6471 0], 'LineStyle', '--');
hold off;
xlim([0,10000])
xlabel('Sample Index');
ylabel('Amplitude');
title('Transmitted Signal with and without AWGN');
legend('show');
grid on;

figure;
% Plot transmitted bitstream without AWGN
subplot(2, 1, 2);
stem(estimatedBitStream(sx:sx+100), 'Marker', 'x', 'DisplayName', 'Estimated Bitstream');
title('Estimated Bitstream');
xlim([0,100])
% Plot estimated bitstream with AWGN
subplot(2, 1, 1);
stem(transmitedBitStream(sx:sx+100), 'Marker', 'o', 'DisplayName', 'Original Bitstream');
title('Original Bitstream');
xlim([0,100])

% Clearing for BER calculations
clearvars
load train;
unquantizedSignal = y;
matchedFilterFlag = 1;
spacing = 100;
noiseVariance = logspace(-3,3,spacing);
Vp = 1;

% Calculating BER of 2PAM
N = 1;
Es = 1;
M = 2^N;
[quantizedSignal,~,~,~,~,~] = MyQuantizer(unquantizedSignal,Vp,N);
transmitedBitStream = MyGraycode(quantizedSignal,Vp,N);

for i = 1:spacing
    receivedSignal = MyMPAM(transmitedBitStream, M, Es);
    receiveSignal = MyAWGNchannel(receivedSignal,noiseVariance(i));
    [~, BER] = DemodulateMPAM(receiveSignal, M, Es, transmitedBitStream, matchedFilterFlag);
    Ber_2_PAM(i) = BER;
end 

% Calculating BER of 2PAM
N = 3;
Es = 1;
M = 2^N;
[quantizedSignal,varLin,varSat,varTeo,SNqR,SNqRTeo] = MyQuantizer(unquantizedSignal,Vp,N);
transmitedBitStream = MyGraycode(quantizedSignal,Vp,N);

for i = 1:spacing
    receivedSignal = MyMPAM(transmitedBitStream, M, Es);
    receiveSignal = MyAWGNchannel(receivedSignal,noiseVariance(i));
    [~, BER] = DemodulateMPAM(receiveSignal, M, Es, transmitedBitStream, matchedFilterFlag);
    Ber_8_PAM(i) = BER;
end 

% Plotting log noise 
figure;
loglog(1:spacing, noiseVariance, 'o-', 'LineWidth', 2);
title('Noise Variance vs. Index of Columns');
xlabel('Index of Columns');
ylabel('Noise Variance');
grid on;

% Plot BER of 2PAM
figure;
semilogx(noiseVariance, Ber_2_PAM, '-o');
title('BER vs Noise Variance for 2-PAM');
xlabel('Noise Variance');
ylabel('Bit Error Rate (BER)');
xlim([0.001,1000]);
ylim([-0.2,0.8]);
grid on;

% Plot BER of 8PAM 
figure;
semilogx(noiseVariance, Ber_8_PAM, '-o');
title('BER vs Noise Variance for 8-PAM');
xlabel('Noise Variance');
ylabel('Bit Error Rate (BER)');
xlim([0.001,1000]);
ylim([-0.2,0.8]);
grid on;

% Functions

% Asignment 1 Functions
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

% Asignment 2 Functions
function transmitSignal = MyMPAM(bitstream,M,Es)
    d = sqrt(3*Es/(M^2-1));
    k = log2(M);
    symbolMatrix = reshape(bitstream, k, length(bitstream)/k)';

    % Map symbols to amplitudes (equally and symmetrically spaced)
    amplitudeLevels = linspace(-(M-1)*d, (M-1)*d, M);
    symbolVector = zeros(1,length(symbolMatrix));
    for i = 1:length(symbolVector)
        for j = 1:k
            symbolVector(i) = symbolVector(i) + symbolMatrix(i,j)*2^-(j-k);
        end
    end 
    transmitSignalfirst = amplitudeLevels(symbolVector+1);
    
    % Rectangular pulse shaping
    sampleRate = 100;
    shapedSignal = rectpulse(transmitSignalfirst, sampleRate);

    % Return shaped signal
    transmitSignal=shapedSignal;
end 

function receiveSignal = MyAWGNchannel(transmitSignal, noiseVariance)
    % Generate white Gaussian noise with specified variance
    noise = wgn(size(transmitSignal, 1), size(transmitSignal, 2), noiseVariance, 'linear');
    
    % Add noise to the transmitted signal
    receiveSignal = transmitSignal + noise;
end

function [estimatedBitStream, BER] = DemodulateMPAM(receivedSignal, M, Es, transmittedBitstream, matchedFilterFlag)
T = 100;
d = sqrt(3*Es/(M^2-1));
k = log2(M);
receivedVector = receivedSignal;
if matchedFilterFlag == 1
    %amplitudeLevels = linspace(-(M-1)*d, (M-1)*d, M);
    segmentSize=floor(length(receivedVector)/T);
    integral=zeros(1,segmentSize);
    for j=1:segmentSize-99
        integral(j)=trapz(receivedSignal((1+((j -1)*T)):(j*T)))/T ;
    end
else
    receivedVector = downsample(receivedSignal, T).';
    integral=receivedVector.';
end

% Estimate symbols from received samples
amplitudeLevels = linspace(-(M-1)*d, (M-1)*d, M);
estimatedSymbols = zeros(1,size(integral, 2));
for i = 1:size(integral, 2)
    [~, index] = min(abs(integral(1, i) - amplitudeLevels));
    estimatedSymbols(i) = index-1;
end

% Convert symbols to bitstream
estimatedSymbolMatrix = de2bi(estimatedSymbols, k, 'left-msb');
estimatedBitStream = reshape(estimatedSymbolMatrix',1, []);

% Calculate Bit Error Rate (BER)
numErrors = sum(estimatedBitStream ~= transmittedBitstream);
BER = numErrors / length(transmittedBitstream);
end

% Asignment 1 Function
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



