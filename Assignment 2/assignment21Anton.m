
M=2^kpv;
Es=1;
k = log2(M);
transmittedSignal = MyMPAM(bitstream, M, Es);





function Signal = MyMPAM(bitstream,M,Es)
    k = log2(M);
    symbolMatrix = reshape(bitstream, k, length(bitstream)/k)';
    symbols=bi2de(symbolMatrix, 'left-msb') + 1;
    % Map symbols to amplitudes (equally and symmetrically spaced)
    amplitudeLevels = linspace(-sqrt(Es), sqrt(Es), M);
    signalAmplitudes = amplitudeLevels(symbols);
    % Create the transmitted signal
    Signal = repelem(signalAmplitudes, k);


end 



function 