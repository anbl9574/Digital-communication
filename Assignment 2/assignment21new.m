clc
clearvars
close all

%Git sug mitt ass

load("Bitstream.mat");
bitstream = estimatedBitStream;
%randi([0 1], 100, 1);
M=4;
Es=4;
k = log2(M);
transmittedSignal = MyMPAM(bitstream, M, Es);
plot(transmittedSignal(1:4000))

function Signal = MyMPAM(bitstream,M,Es)
    k = log2(M);
    symbolMatrix = reshape(bitstream, k, length(bitstream)/k)';
    symbols=bi2de(symbolMatrix, 'left-msb') + 1;
    % Map symbols to amplitudes (equally and symmetrically spaced)
    amplitudeLevels = linspace(-sqrt(Es), sqrt(Es), M);
    %disp(amplitudeLevels)
    signalAmplitudes = amplitudeLevels(symbols);
    %disp(signalAmplitudes)
    % Create the transmitted signal
    Signal = repelem(signalAmplitudes, k);
%     for i = 1:length(Signal)
%         Es_calc =+ Signal(i)^2;
%     end
%     disp(Es_calc)
end 