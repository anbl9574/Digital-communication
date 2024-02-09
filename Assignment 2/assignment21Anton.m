clc
clearvars
close all

load("Bitstream.mat");
bitstream = estimatedBitStream;

transmittedBitstream=bitstream;
%randi([0 1], 100, 1);
M=4;
Es=1;
T=100;
receivedSignal = MyMPAM(bitstream, M, Es);
plot(receivedSignal(1:10000))%transmitted
[estimatedBitstream, BER] = DemodulateMPAM(receivedSignal,M,Es,transmittedBitstream)

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
         estimatedSymbols = zeros(size(receivedMatrix, 1), 1);
         for i = 1:size(receivedMatrix, 1)
            [~, index] = min(receivedMatrix(i, 1:4) -  linspace(-(M-1)*d, (M-1)*d, M)); %Something is wrong here
            estimatedSymbols(i) = index - 1;  

            % Convert symbols to bitstream
            estimatedSymbolMatrix = de2bi(estimatedSymbols, k, 'left-msb');
            estimatedBitstream = reshape(estimatedSymbolMatrix', 1, []);

            % Calculate Bit Error Rate (BER)
            %numErrors = sum(estimatedBitstream ~= transmittedBitstream); %Wrong here
            %BER = numErrors / length(transmittedBitstream);

         end
end 
