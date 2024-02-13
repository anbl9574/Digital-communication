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
receivedSignal = MyMPAM(bitstream, M, Es);
noiseVariance = logspace (0,3,size(receivedSignal,2)) ;
receiveSignal = MyAWGNchannel(receivedSignal,noiseVariance);
[estimatedBitstream, BER] = DemodulateMPAM(receivedSignal,M,Es,transmittedBitstream);


%%
%Plots

plot(receivedSignal(1:10000))%transmitted

plot(receiveSignal(1:10000))
% Plotting log noise 
figure;
loglog(1:size(receivedSignal, 2), noiseVariance, 'o-', 'LineWidth', 2);
title('Noise Variance vs. Index of Columns');
xlabel('Index of Columns');
ylabel('Noise Variance');
grid on;

% Generate a random starting index
startIndex = randi(length(bitstream) -100 + 1);
% Select a range of 50 elements
endIndex = startIndex + 99;
subsetBitstream = bitstream(startIndex:endIndex);
subsetEstimatedBitstream = estimatedBitstream(startIndex:endIndex);
figure;

subplot(2, 1, 1);
stem(subsetBitstream, 'Marker', 'o', 'DisplayName', 'Original Bitstream');
title('Original Bitstream');

subplot(2, 1, 2);
stem(subsetEstimatedBitstream, 'Marker', 'x', 'DisplayName', 'Estimated Bitstream');
title('Estimated Bitstream');
legend('show');


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
    %Signal = repelem(amplitudeLevels, k);
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

         receivedMatrix = reshape(receivedSignal, T, length(receivedSignal)/T)';
         amplitudeLevels = linspace(-(M-1)*d, (M-1)*d, M);

       % Estimate symbols from received samples
        [~, symbols] = ismember(receivedMatrix, amplitudeLevels);
        

        symbols=symbols(:,1)-Es;
        symbols = symbols(symbols ~= -1);

        % Convert symbols to bitstream
        estimatedSymbolMatrix = de2bi(symbols, k, 'left-msb');
        estimatedBitstream = reshape(estimatedSymbolMatrix', 1, []);

        % Calculate Bit Error Rate (BER)
        numErrors = sum(estimatedBitstream ~= transmittedBitstream);
        BER = numErrors / length(transmittedBitstream);

end


function receiveSignal = MyAWGNchannel(transmitSignal, noiseVariance)
    % Add white Gaussian noise to the transmitted signal
    noise = sqrt(noiseVariance) .* randn(size(transmitSignal));

    receiveSignal = transmitSignal + noise;
end



