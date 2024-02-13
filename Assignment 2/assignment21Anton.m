clc
clearvars
close all

%%
%Defining parameters and loading in bitstream from assignment 2

load("Bitstream.mat");
bitstream = estimatedBitStream;
transmittedBitstream=bitstream;
M=4;
Es=1;
T=100;
%%

%Calling the functions
spacing=100;

Ber_2_PAM=zeros(spacing,1);
Ber_4_PAM=zeros(spacing,1);

receivedSignal = MyMPAM(bitstream, M, Es);
%noiseVariance = logspace (0,3) ;
receiveSignal = MyAWGNchannel(receivedSignal,10);

[estimatedBitstream, BER] = DemodulateMPAM(receiveSignal,M,Es,transmittedBitstream);
noiseVariance = logspace(-1,-3,spacing);

%testing
for i = 1:spacing
    receivedSignal = MyMPAM(bitstream, M, Es);
    receiveSignal = MyAWGNchannel(receivedSignal,noiseVariance(i));
    [estimatedBitstream, BER] = DemodulateMPAM(receiveSignal,M,Es,transmittedBitstream);
    % Store BER for plotting
    Ber_2_PAM(i) = BER;
end 

M=4;
Es=1;
T=100;


for i = 1:spacing
    receivedSignal = MyMPAM(bitstream, M, Es);
    receiveSignal = MyAWGNchannel(receivedSignal,noiseVariance(i));
    [estimatedBitstream, BER] = DemodulateMPAM(receiveSignal,M,Es,transmittedBitstream);
    % Store BER for plotting
    Ber_4_PAM(i) = BER;
end 



%%
%Plots

plot(receivedSignal(1:10000))%transmitted

plot(receiveSignal(1:10000))
% Plotting log noise 
%figure;
%loglog(1:size(receivedSignal, 2), noiseVariance, 'o-', 'LineWidth', 2);
%title('Noise Variance vs. Index of Columns');
%xlabel('Index of Columns');
%ylabel('Noise Variance');
%grid on;

% Generate a random starting index
startIndex = randi(length(bitstream) -100 + 1);
% Select a range of 50 elements
endIndex = startIndex + 99;
subsetBitstream = bitstream(startIndex:endIndex);
subsetEstimatedBitstream = estimatedBitstream(startIndex:endIndex);
figure;


%Plot for estimated signal vs original bitstream
subplot(2, 1, 1);
stem(subsetBitstream, 'Marker', 'o', 'DisplayName', 'Original Bitstream');
title('Original Bitstream');

subplot(2, 1, 2);
stem(subsetEstimatedBitstream, 'Marker', 'x', 'DisplayName', 'Estimated Bitstream');
title('Estimated Bitstream');
legend('show');

% Plot
figure;
semilogx(noiseVariance, Ber_2_PAM, '-o');
title('BER vs Noise Variance for 2-PAM');
xlabel('Noise Variance');
ylabel('Bit Error Rate (BER)');
grid on;

% Plot
figure;
semilogx(noiseVariance, Ber_4_PAM, '-o');
title('BER vs Noise Variance for 4-PAM');
xlabel('Noise Variance');
ylabel('Bit Error Rate (BER)');
grid on;

%%
%Functions

function receivedSignal = MyMPAM(bitstream,M,Es)
    d = sqrt(3*Es/(M^2-1));
    k = log2(M);
    symbolMatrix = reshape(bitstream, k, length(bitstream)/k)';
    symbols=bi2de(symbolMatrix, 'left-msb') + Es;
    % Map symbols to amplitudes (equally and symmetrically spaced)
    amplitudeLevels = linspace(-(M-1)*d, (M-1)*d, M);
    %disp(amplitudeLevels)
    Signaltransmit = amplitudeLevels(symbols);
    %disp(signalAmplitudes)
    % Create the transmitted signal
    E = 0;
    for i = 1:length(amplitudeLevels)
        E = E + amplitudeLevels(i)^2;
    end
    E = E/4;
    if length(Signaltransmit)==length(bitstream)
        transmitSignal=Signaltransmit;
    else 
        transmitSignal=padarray(Signaltransmit,[0,length(bitstream)-length(Signaltransmit)],'post');
    end 
    % Rectangular pulse shaping
    T=100;
    sampleRate = T;
    shapedSignal = rectpulse(transmitSignal, sampleRate);

    % Return shaped signal
    receivedSignal=shapedSignal;


end 

function [estimatedBitstream, BER] = DemodulateMPAM(receivedSignal, M, Es, transmittedBitstream, T)

        % Reshape received signal into matrix with T elements per rowT=100;
         T=100;
         d = sqrt(3*Es/(M^2-1));
         k = log2(M);

         receivedMatrix = downsample(receivedSignal, T).';

         %receivedMatrix = reshape(receivedSignal, T, length(receivedSignal)/T)';
         amplitudeLevels = linspace(-(M-1)*d, (M-1)*d, M);

       % Estimate symbols from received samples
         estimatedSymbols = zeros(size(receivedMatrix, 1), 1);

        for i = 1:size(receivedMatrix, 1)
        [~, index] = min(abs(receivedMatrix(i, 1) - amplitudeLevels));
        estimatedSymbols(i) = index-1;  % Adjust for 0-based indexing
        end
        
       % symbols = estimatedSymbols(estimatedSymbols ~= -1);

        % Convert symbols to bitstream
        estimatedSymbolMatrix = de2bi(estimatedSymbols, k, 'left-msb');
        estimatedBitstream = reshape(estimatedSymbolMatrix', k, []);
        
        % Calculate Bit Error Rate (BER)

        disp(size(estimatedBitstream))
        disp(size(transmittedBitstream))

        numErrors = sum(estimatedBitstream ~= transmittedBitstream);
        BER = numErrors / length(transmittedBitstream);
end


function receiveSignal = MyAWGNchannel(transmitSignal, noiseVariance)
    % Add white Gaussian noise to the transmitted signal
    noise = sqrt(noiseVariance) * randn(size(transmitSignal));

    receiveSignal = transmitSignal + noise;
end



