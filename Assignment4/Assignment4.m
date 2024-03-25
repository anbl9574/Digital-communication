close all
clearvars
clc

load("Bitstream4bit.mat");
c = estimatedBitStream;

HammingTest(c)
BERplot(c)

%Functions
% function c = OwnEncodeHAM74(b)
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
    [~,G] = hammgen(3);
    for i = 1:length(b)/4
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
    [~,G] = hammgen(4);
    for i = 1:length(b)/11
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
    H = hammgen(3);
    for i = 1:length(c)/7
        d = c(7*(i-1)+1:7*(i-1)+7);
        s = d*H';
        for j = 1:3
            s(j) = mod(s(j),2);
        end
        if isequal(s,[0,0,0])
            bEst = [bEst,d(4:7)];
        else
            [~,idx] = ismember(s,H','rows');
            d(idx) = ~d(idx);
            bEst = [bEst,d(4:7)];
        end
    end
end

function bEst = decodeHAM1511(c)
    % Decoding
    bEst = [];
    H = hammgen(4);
    for i = 1:length(c)/15
        d = c(15*(i-1)+1:15*(i-1)+15);
        s = d*H';
        for j = 1:4
            s(j) = mod(s(j),2);
        end
        if isequal(s,[0,0,0,0])
            bEst = [bEst,d(5:15)];
        else
            [~,idx] = ismember(s,H','rows');
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

function HammingTest(c)
    c1 = encodeHAM74(c);
    c2 = encodeHAM1511(c);
    e1 = decodeHAM74(c1);
    e2 = decodeHAM1511(c2);
    e2 = e2(1:length(c));
    
    x4 = (1:44);
    x7 = (1:77);
    x11 = (1:44);
    x15 = (1:60);
    
    figure;
    % First subplot
    subplot(3, 1, 1);
    stairs(x4, c(1:44), 'b');
    title('Original Bits');
    xlabel('Bits');
    ylabel('Value');
    ylim([-0.1,1.1])
    xlim([1,44])

    % Second subplot
    subplot(3, 1, 2);
    stairs(x7, c1(1:77), 'r');
    title('Hamming (7,4) Coded Bits');
    xlabel('Bits');
    ylabel('Value');
    ylim([-0.1,1.1])
    xlim([1,77])

    % Third subplot
    subplot(3, 1, 3);
    stairs(x4, e1(1:44), 'b');
    title('Decoded Bits');
    xlabel('Bits');
    ylabel('Value');
    ylim([-0.1,1.1])
    xlim([1,44])
    
    figure;
    title('Hamming(15,11)');
    % First subplot
    subplot(3, 1, 1);
    stairs(x11, c(1:44), 'b');
    title('Original Bits');
    xlabel('Bits');
    ylabel('Value');
    ylim([-0.1,1.1])
    xlim([1,44])

    % Second subplot
    subplot(3, 1, 2);
    stairs(x15, c2(1:60), 'r');
    title('Hamming (15,11) Coded Bits');
    xlabel('Bits');
    ylabel('Value');
    ylim([-0.1,1.1])
    xlim([1,60])

    % Third subplot
    subplot(3, 1, 3);
    stairs(x11, e2(1:44), 'b');
    title('Decoded Bits');
    xlabel('Bits');
    ylabel('Value');
    ylim([-0.1,1.1])
    xlim([1,44])
end

function BERplot(c)
%c = c(1:1000);
%c = [c,c,c,c];
%c = [c,c,c,c];

%Theoretical BPSK BER
Eb_No = logspace(-1,3,100);
Eb = 1;
No = Eb./Eb_No;
NO = 10.^(No / 10);
BERBPSKTHEO = qfunc(sqrt(2*Eb_No));


%Theoretical BPSK with HAM74 BER
n = 7;
k = 4;
p = qfunc(sqrt(2*Eb_No));
BERBPSKTHEO7 = p-p.*(1-p).^(n-1);
%0.5.*erfc(sqrt((k/n).*Eb_No)).*(1- (1-0.5.*erfc(sqrt((k/n).*Eb_No))).^(n-1));

%Theoretical BPSK with HAM1511 BER
n = 15;
k = 11;
BERBPSKTHEO15 = p-p.*(1-p).^(n-1);
%qfunc(sqrt(2*Eb_No)).*(1 - (1-0.5.*erfc(sqrt((k/n).*Eb_No))).^(n-1));

%Empirical BPSK with HAM74 BER
BERBPSK7 = (1:length(No));
e = encodeHAM74(c);
s = mapontoBPSK(e);
for i = 1:length(No)     
    n = channel2(s,No(i));
    h = detectBPSK(n);
    b = decodeHAM74(h);
    % Claculate and store BER for plotting
    BER = sum(abs(c - b))/length(c);
    BERBPSK7(i) = BER;
end

BERBPSK15 = (1:length(No));
e = encodeHAM1511(c);
s = mapontoBPSK(e);
for i = 1:length(No)
    r = channel2(s,No(i));
    h = detectBPSK(r);
    b = decodeHAM1511(h);
    b = b(1:length(c));
    
    % Claculate and store BER for plotting
    BER = sum(c ~= b)/length(c);
    BERBPSK15(i) = BER;
end

% Plotting Theoretical BER
figure;
semilogy(db(Eb_No,'power'),BERBPSKTHEO7);
hold on;
semilogy(db(Eb_No,'power'),BERBPSKTHEO15);
semilogy(db(Eb_No,'power'),BERBPSKTHEO);
title('Theoretical BER for BPSK using HAM74 and HAM1511')
xlabel('Signal Noise Ratio (SNR [dB])');
ylabel('Bit Error Rate (BER)');
ylim([10^-7,1]);
grid on;
legend('Hamming (7,4)','Hamming (15,11)', 'Uncoded BPSK', 'Location', 'southwest');

% Plotting Empirical BER for HAM74
figure;
semilogy(db(Eb_No,'power'),BERBPSK7, 'bx');
hold on;
semilogy(db(Eb_No,'power'),BERBPSKTHEO7);
semilogy(db(Eb_No,'power'),BERBPSKTHEO);
title('Empirical BER for BPSK using HAM74');
xlabel('Signal Noise Ratio (SNR [dB])');
ylabel('Bit Error Rate (BER)');
ylim([10^-7,1]);
grid on;
legend('Empirical Hamming(7,4)', 'Theoretical Hamming(7,4)', 'Uncoded BPSK', 'Location', 'southwest');

% Plotting Empirical BER for HAM1511
figure;
semilogy(db(Eb_No,'power'),BERBPSK15, 'bx');
hold on;
semilogy(db(Eb_No,'power'),BERBPSKTHEO15);
semilogy(db(Eb_No,'power'),BERBPSKTHEO);
title('Empirical BER for BPSK using HAM1511');
xlabel('Signal Noise Ratio (SNR [dB])');
ylabel('Bit Error Rate (BER)');
ylim([10^-7,1]);
grid on;
legend('Empirical Hamming(15,11)', 'Theoretical Hamming(15,11)', 'Uncoded BPSK', 'Location', 'southwest');


end
