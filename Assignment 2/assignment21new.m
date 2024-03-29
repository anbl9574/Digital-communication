clc
clearvars
close all

load("Bitstream.mat");
bitstream = estimatedBitStream;
%randi([0 1], 100, 1);
M=4;
Es=1;
%transmittedSignal = MyMPAM(bitstream, M, Es);

%function Signal = MyMPAM(bitstream,M,Es)
    d = sqrt(3*Es/(M^2-1))
    k = log2(M);
    symbolMatrix = reshape(bitstream, k, length(bitstream)/k)';
    symbols=bi2de(symbolMatrix, 'left-msb') + 1;
    % Map symbols to amplitudes (equally and symmetrically spaced)
    amplitudeLevels = linspace(-(M-1)*d, (M-1)*d, M);
    %disp(amplitudeLevels)
    Signal = amplitudeLevels(symbols);
    %disp(signalAmplitudes)
    % Create the transmitted signal
    %Signal = repelem(amplitudeLevels, k);
    E = 0;
    for i = 1:length(amplitudeLevels)
        E = E + amplitudeLevels(i)^2;
    end
    E = E/4
%end 
plot(Signal(1:500))%transmitted


