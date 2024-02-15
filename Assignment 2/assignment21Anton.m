clc
clearvars
close all

%%
%Defining parameters and loading in bitstream from assignment 2

load("Bitstream1bit.mat");
bitstream = estimatedBitStream;

M = 2;
Es = 1;
T=100;
matchedFilterFlag=1;
%%

%Calling the functions
spacing=100;

Ber_2_PAM=zeros(spacing,1);
Ber_8_PAM=zeros(spacing,1);

transmitSignal = MyMPAM(bitstream, M, Es);
receiveSignal = MyAWGNchannel(transmitSignal,0.01);
%%

% Plotting the signal with noise against the one without
figure;

% Plot received signal with AWGN
plot(receiveSignal(1:6700), 'LineWidth', 2, 'DisplayName', 'Received Signal with AWGN');
hold on;

% Plot transmitted signal without AWGN with transparency
plot(transmitSignal(1:6700), 'LineWidth', 2, 'DisplayName', 'Transmitted Signal', 'Color',[1 0.6471 0], 'LineStyle', '--');

hold off;

xlabel('Sample Index');
ylabel('Amplitude');
title('Transmitted Signal with and without AWGN');
legend('show');

grid on;
% Generate a random starting index
startIndex = randi(length(bitstream) -100 + 1);
% Select a range of 50 elements
endIndex = startIndex + 99;
subsetBitStream = bitstream(startIndex:endIndex);
subsetEstimatedBitstream = estimatedBitStream(startIndex:endIndex);
figure;


%Plot for estimated signal vs original bitstream
subplot(2, 1, 1);
stem(subsetBitStream, 'Marker', 'o', 'DisplayName', 'Original Bitstream');
title('Original Bitstream');

subplot(2, 1, 2);
stem(subsetEstimatedBitstream, 'Marker', 'x', 'DisplayName', 'Estimated Bitstream');
title('Estimated Bitstream');
legend('show');
%%

transmittedBitstream=bitstream;
[estimatedBitStream, BER] = DemodulateMPAM(receiveSignal, M, Es, transmittedBitstream, matchedFilterFlag);
noiseVariance = logspace(-3,-1,spacing);

clear receiveSignal;

%testing
for i = 1:spacing
    receivedSignal = MyMPAM(transmittedBitstream, M, Es);
    receiveSignal = MyAWGNchannel(receivedSignal,noiseVariance(i));
    [estimatedBitStream, BER] = DemodulateMPAM(receiveSignal, M, Es, transmittedBitstream, matchedFilterFlag);
    % Store BER for plotting
    Ber_2_PAM(i) = BER;
end 

clear bitstream;
clear estimatedBitStream;
clear receivedSignal;
clear receiveSignal;
clear transmittedBitstream;
load("Bitstream3bit.mat");
transmittedBitstream=estimatedBitStream;

M=8;
Es=1;
T=100;

for i = 1:spacing
    receivedSignal = MyMPAM(transmittedBitstream, M, Es);
    receiveSignal = MyAWGNchannel(receivedSignal,noiseVariance(i));
    [estimatedBitStream, BER] = DemodulateMPAM(receiveSignal, M, Es, transmittedBitstream, matchedFilterFlag);
    % Store BER for plotting
    Ber_8_PAM(i) = BER;
end 




%%
%Plots

%plot(receiveSignal(1:500000))%with AWGN
% Plotting log noise 
figure;
loglog(1:spacing, noiseVariance, 'o-', 'LineWidth', 2);
title('Noise Variance vs. Index of Columns');
xlabel('Index of Columns');
ylabel('Noise Variance');
grid on;


% Plot
figure;
semilogx(noiseVariance, Ber_2_PAM, '-o');
title('BER vs Noise Variance for 2-PAM');
xlabel('Noise Variance');
ylabel('Bit Error Rate (BER)');
grid on;

% Plot
figure;
semilogx(noiseVariance, Ber_8_PAM, '-o');
title('BER vs Noise Variance for 8-PAM');
xlabel('Noise Variance');
ylabel('Bit Error Rate (BER)');
grid on;

%%
%Functions

function transmitSignal = MyMPAM(bitstream,M,Es)
    d = sqrt(3*Es/(M^2-1));
    
    k = log2(M);

    if mod(length(bitstream), k) ~= 0

    % Calculate the number of zeros needed for padding
    paddingSize = k - mod(length(bitstream), k);
    % Pad the bitstream with zeros
    bitstream = transpose([bitstream.', zeros(1, paddingSize)]);   
    else 
        bitstream=bitstream;
    end 
    symbolMatrix = reshape(bitstream, k, length(bitstream)/k)';
    
    % Map symbols to amplitudes (equally and symmetrically spaced)

    amplitudeLevels = linspace(-(M-1)*d, (M-1)*d, M);
    if k>1
    symbols=bi2de(symbolMatrix,k-1,'left-msb') +Es;
    else 
        symbols=symbolMatrix+Es;

    end 
    % Map symbols to amplitudes (equally and symmetrically spaced)
    %disp(amplitudeLevels)
    symbols(symbols > M) = M; %Don't know
    transmitSignalfirst = amplitudeLevels(symbols);
    %disp(signalAmplitudes)
    % Create the transmitted signal
    E = 0;
    for i = 1:length(amplitudeLevels)
        E = E + amplitudeLevels(i)^2;
    end
    E = E/4;

    % Rectangular pulse shaping
    T=100;
    sampleRate = T;
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

T=100;
d = sqrt(3*Es/(M^2-1));
k = log2(M);
receivedMatrix=receivedSignal;

if matchedFilterFlag==1
    amplitudeLevels = linspace(-(M-1)*d, (M-1)*d, M);

    segmentSize=floor(length(receivedMatrix)/T);
    integral1=zeros(1,segmentSize);
    for j=1:segmentSize-99
        integral1(j)=trapz(receivedSignal((1+((j -1)*T)):(j*T)))/T ;
    end
else
    receivedMatrix = receivedMatrix;
    receivedMatrix = downsample(receivedSignal, T).';
    integral1=receivedMatrix.';

end
amplitudeLevels = linspace(-(M-1)*d, (M-1)*d, M);
% Estimate symbols from received samples

estimatedSymbols = zeros(1,size(integral1, 2));

for i = 1:size(integral1, 2)
    [~, index] = min(abs(integral1(1, i) - amplitudeLevels));
    estimatedSymbols(i) = index-1;
end

% Convert symbols to bitstream
estimatedSymbolMatrix = de2bi(estimatedSymbols, k, 'left-msb');
estimatedBitStream = reshape(estimatedSymbolMatrix',1, []);
estimatedBitStream=estimatedBitStream(1:size(transmittedBitstream,2));
% Calculate Bit Error Rate (BER)

numErrors = sum(estimatedBitStream ~= transmittedBitstream);
BER = numErrors / length(transmittedBitstream);
end



