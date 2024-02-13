clc
clearvars
close all

%%
%Defining parameters and loading in bitstream from assignment 2

load("Bitstream1bit.mat");
bitstream = estimatedBitStream;
transmittedBitstream=bitstream;
M=2;
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
noiseVariance = logspace(-3,-1,spacing);

%testing
for i = 1:spacing
    receivedSignal = MyMPAM(bitstream, M, Es);
    receiveSignal = MyAWGNchannel(receivedSignal,noiseVariance(i));
    [estimatedBitstream, BER] = DemodulateMPAM(receiveSignal,M,Es,transmittedBitstream);
    % Store BER for plotting
    Ber_2_PAM(i) = BER;
end 



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

    if mod(length(bitstream), k) ~= 0

    % Calculate the number of zeros needed for padding
    paddingSize = k - mod(length(bitstream), k);
    % Pad the bitstream with zeros
    bitstream = padarray(bitstream, [0, paddingSize], 'post');   
    else 
        bitstream=bitstream;
    end 
    symbolMatrix = reshape(bitstream, k, length(bitstream)/k)';
    
    % Map symbols to amplitudes (equally and symmetrically spaced)

    amplitudeLevels = linspace(-(M-1)*d, (M-1)*d, M);
    if k>1
    symbols=bi2de(symbolMatrix,k,'left-msb') +Es;
    else 
        symbols=symbolMatrix+Es;

    end 
    % Map symbols to amplitudes (equally and symmetrically spaced)
    %disp(amplitudeLevels)
    symbols(symbols > M) = M; %Don't know
    transmitSignal = amplitudeLevels(symbols);
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
        estimatedSymbols(i) = index-1; 
        end
        
        % Convert symbols to bitstream
        estimatedSymbolMatrix = de2bi(estimatedSymbols, k, 'left-msb');
        estimatedBitstream = reshape(estimatedSymbolMatrix', 1, []);
        estimatedBitstream=estimatedBitstream(1:size(transmittedBitstream,2));
        % Calculate Bit Error Rate (BER)
    
        numErrors = sum(estimatedBitstream ~= transmittedBitstream);
        BER = numErrors / length(transmittedBitstream);
end


function receiveSignal = MyAWGNchannel(transmitSignal, noiseVariance)
    % Add white Gaussian noise to the transmitted signal
    noise = sqrt(noiseVariance) * randn(size(transmitSignal));

    receiveSignal = transmitSignal + noise;
end



