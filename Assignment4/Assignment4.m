close all
clearvars
clc

load("Bitstream4bit.mat");
c = estimatedBitStream;

c1 = encodeHAM74(c);
c2 = encodeHAM1511(c);
e1 = decodeHAM74(c1);
e2 = decodeHAM1511(c2);
e2 = e2(1:length(c));
est = [c;e1;e2];

ie1 = isequal(c,e1);
ie2 = isequal(c,e2);

BERplot(c)

%Functions
% function c = OWNencodeHAM74(b)
%     % Padding
%     p = mod(length(b),4);
%     if (p > 0)
%         z = zeros(1,4-p);
%         b = [b,z];
%     end
%     % Encoding
%     c = [];
%     for i = 1:length(b)/4
%         t = zeros(1,7);
%         t(4) = b(4*(i-1)+1);
%         t(5) = b(4*(i-1)+2);
%         t(6) = b(4*(i-1)+3);
%         t(7) = b(4*(i-1)+4);
%         t(1) = xor(xor(t(4),t(6)),t(7));
%         t(2) = xor(xor(t(4),t(5)),t(6));
%         t(3) = xor(xor(t(5),t(6)),t(7));
%         c = [c,t];
%     end
% end

function c = encodeHAM74(b)
    % Padding
    p = mod(length(b),4);
    if (p > 0)
        z = zeros(1,4-p);
        b = [b,z];
    end
    % Encoding
    c = [];
    for i = 1:length(b)/4
        [~,G] = hammgen(3);
        d = b(4*(i-1)+1:4*(i-1)+4);
        t = d*G;
        for j = 1:7 
            t(j) = mod(t(j),2);
        end
        c = [c,t];
    end
end

function c = encodeHAM1511(b)
    % Padding
    p = mod(length(b),11);
    if (p > 0)
        z = zeros(1,11-p);
        b = [b,z];
    end
    % Encoding
    c = [];
    for i = 1:length(b)/11
        [~,G] = hammgen(4);
        d = b(11*(i-1)+1:11*(i-1)+11);
        t = d*G;
        for j = 1:15
            t(j) = mod(t(j),2);
        end
        c = [c,t];
    end
end

function bEst = decodeHAM74(c)
    % Decoding
    bEst = [];
    for i = 1:length(c)/7
        H = hammgen(3);
        d = c(7*(i-1)+1:7*(i-1)+7);
        s = d*H';
        for j = 1:3
            s(j) = mod(s(j),2);
        end
        if isequal(s,[0,0,0])
            bEst = [bEst,d(4:7)];
        else
            [~,idx] = ismember(s,H');
            d(idx) = ~d(idx);
            bEst = [bEst,d(4:7)];
        end
    end
end

function bEst = decodeHAM1511(c)
    % Decoding
    bEst = [];
    for i = 1:length(c)/15
        H = hammgen(4);
        d = c(15*(i-1)+1:15*(i-1)+15);
        s = d*H';
        for j = 1:4
            s(j) = mod(s(j),2);
        end
        if isequal(s,[0,0,0,0])
            bEst = [bEst,d(5:15)];
        else
            [~,idx] = ismember(s,H');
            d(idx) = ~d(idx);
            bEst = [bEst,d(5:15)];
        end
    end
end

function noisySignal = channel2(transmittedSignal,noisePower)
% Calculate the variance
variance = noisePower/2;

% Noise for real and imaginary parts
noise_real = sqrt(variance) * randn(size(transmittedSignal));
noise_imaginary = sqrt(variance) * randn(size(transmittedSignal));

% Combine to form complex noise
noise = noise_real + 1i * noise_imaginary;

% Add noise to transmitted signal
noisySignal = transmittedSignal + noise;
end

function s = mapontoBPSK(c)
% Convert bits to symbols, BPSK mapping: 0 -> -1, 1 -> 1
s = 2*c - 1; 
end

function c = detectBPSK(r)
% Decision threshold for BPSK demodulation
threshold = 0;

% Perform detection
c = (1:length(r));
for i = 1:length(r)
    if (r(i) > threshold)
        c(i) = 1;
    else
        c(i) = 0;
    end
end
end


function BERplot(c)
% Repeating vector for larger amount of samples
%c = [c,c,c,c];%,c,c,c,c,c,c,c,c,c,c,c,c];

% Noise variance vector
No = logspace(-1.5,1,15);

%testing

BERBPSK7 = (1:length(No));
for i = 1:length(No)
    e = encodeHAM74(c);
    s = mapontoBPSK(e);
    r = channel2(s,No(i));
    h = detectBPSK(r);
    b = decodeHAM74(h);
    
    % Claculate and store BER for plotting
    BER = sum(c ~= b)/length(c);
    BERBPSK7(i) = BER;
end

BERBPSK15 = (1:length(No));
for i = 1:length(No)
    e = encodeHAM1511(c);
    s = mapontoBPSK(e);
    r = channel2(s,No(i));
    h = detectBPSK(r);
    b = decodeHAM1511(h);
    b = b(1:length(c));
    
    % Claculate and store BER for plotting
    BER = sum(c ~= b)/length(c);
    BERBPSK15(i) = BER;
end

% Calculating theoretical BER
Es = 1;
Eb = Es/log2(16);
EsNo = Es./No;
EbNo = Eb./No;
EbNo_dB = 10*log10(EbNo);

p = qfunc(sqrt(2*EbNo));

BERBPSKTHEO7 = p-p.*(1-p).^2;%qfunc(sqrt(2*EsNo));
BERBPSKTHEO15 = p-p.*(1-p).^3;

% Plotting BER
figure;
semilogy(EbNo_dB, BERBPSK7, '.-');
hold on;
semilogy(EbNo_dB, BERBPSK15, '.-');
semilogy(EbNo_dB, BERBPSKTHEO7, '.-');
semilogy(EbNo_dB, BERBPSKTHEO15, '.-');
title('BER vs SNR for BPSK');
xlabel('Signal Noise Ratio (SNR [dB])');
ylabel('Bit Error Rate (BER)');
ylim([10^-5,1]);
grid on;
legend('Hamming(7,4)','Hamming(15,11)', 'Theoretical BPSK', 'Location', 'southwest');
end
