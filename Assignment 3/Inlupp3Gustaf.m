close all
clearvars
clc

load("Bitstream4bit.mat");
c = estimatedBitStream;

TEST()
QPSK()
TEST_16QAM(c)
BERCURVES(c)

%Functions

function noisySignal = channel2(transmittedSignal,noisePower)
% Calculate the variance of the real and imaginary parts of noise
variance = noisePower / 2;

% Generate white Gaussian noise for real and imaginary parts
noise_real = sqrt(variance) * randn(size(transmittedSignal));
noise_imaginary = sqrt(variance) * randn(size(transmittedSignal));

% Combine real and imaginary parts to form complex noise
noise = noise_real + 1i * noise_imaginary;

% Add noise to transmitted signal
noisySignal = transmittedSignal + noise;
end

function QPSK()
% Number of symbols
numSymbols = 4;

% Specify phase angles for QPSK symbols
phase_angles = [pi/4, 3*pi/4, 5*pi/4, 7*pi/4]; % Phase angles for QPSK symbols

% Create complex symbols with unit magnitude and specified phases
symbols = (1:4);
for i = 1:numSymbols
    symbols(i) = exp(1i * phase_angles(i));
end

% Verify average symbol power
average_power = mean(abs(symbols).^2);
disp(['QPSK Average symbol power: ', num2str(average_power)]);

% Plot the symbols and detection regions in the complex plane
figure;
plot(real(symbols), imag(symbols), 'cx', 'MarkerSize', 10,'LineWidth', 3);
xlabel('In-phase');
ylabel('Quadrature');
title('QPSK Symbols with Detection Regions');

xlim([-3.2, 3.2]);
ylim([-3.2, 3.2]);
axis square;
grid on;
hold on;

% Define detection regions and plot
xline(0, 'r--', 'LineWidth', 3);
yline(0, 'r--', 'LineWidth', 3);
legend('Sybols', 'Decition Boundries', 'Location', 'southeast')
    
hold off;

% Verify for different noise power
SignalSymbols = symbols(2) * ones(1, 1000);
noisePower2 = [0.01, 0.1, 1];

for i = 1:length(noisePower2)
    ReccivedSymbols = channel2(SignalSymbols,noisePower2(i));
    figure;
    plot(real(symbols), imag(symbols), 'cx', 'MarkerSize', 10, 'LineWidth', 3);
    hold on;
    plot(real(ReccivedSymbols), imag(ReccivedSymbols), 'b.', 'MarkerSize', 10, 'LineWidth', 2);
    plot(real(SignalSymbols), imag(SignalSymbols), 'rx', 'MarkerSize', 10, 'LineWidth', 3);
    xlabel('Real-Part');
    ylabel('Imaginary-Part');
    
    title(['QPSK Sent Vs Reccived Symbols with Detection Regions, noise variance: ', num2str(noisePower2(i))]);
    xline(0, 'r--', 'LineWidth', 3);
    yline(0, 'r--', 'LineWidth', 3);
    
%     xmax = max(real(ReccivedSymbols));
%     xmin = min(real(ReccivedSymbols));
%     xlims = max(abs(xmax),abs(xmin));
%     lim = xlims + 0.1;
    xlim([-3.2, 3.2]);
    ylim([-3.2, 3.2]);
    legend('Sybols', 'Noisy Reccived Symbols', 'Sent Symbol', 'Decition Boundries', 'Location', 'southeast')
    axis square;
    grid on;
    hold off;
    
    %SER
    ErrorSymbol = sum(angle(ReccivedSymbols) > pi | angle(ReccivedSymbols) < pi/2);
    SER = ErrorSymbol/length(ReccivedSymbols);
    disp(['QPSK Noise varriance: (', num2str(noisePower2(i)),'), Empirical SER: ', num2str(SER)]);%, ', Empirical BER: ' num2str(SER)]);
    
end
end

function s = maponto16QAM(c)
Es = 1;
% Create 16QAM Symbols
P = 4*Es/(sqrt(2)+2*sqrt(10)+sqrt(18));
Symbols = [-3+3i,-3+1i,-3-3i,-3-1i, ...
           -1+3i,-1+1i,-1-3i,-1-1i, ...
            3+3i, 3+1i, 3-3i, 3-1i, ...
            1+3i, 1+1i, 1-3i, 1-1i]*P;
        
% Map bits to 16QAM symbols
s = (1:length(c)/4);
for i = 1:length(c)/4
    a = c(4*(i-1)+1:4*(i-1)+4);
    s(i) = Symbols(bin2dec(char(a + '0'))+1);
end
end

function s = mapontoBPSK(c)
% Convert bits to symbols
s = 2*c - 1; % BPSK mapping: 0 -> -1, 1 -> 1
end

function c = detect16QAM(r)
Es = 1;
% Define 16-QAM constellation points
P = 4*Es/(sqrt(2)+2*sqrt(10)+sqrt(18));
Symbols = [-3+3i,-3+1i,-3-3i,-3-1i, ...
           -1+3i,-1+1i,-1-3i,-1-1i, ...
            3+3i, 3+1i, 3-3i, 3-1i, ...
            1+3i, 1+1i, 1-3i, 1-1i]*P;

% Perform detection
c = [];
for i = 1:length(r)
    % Compute Euclidean distances
    distances = abs(r(i) - Symbols).^2;
    [~, index] = min(distances); % Find index of the minimum distance
    
    % Map index to binary sequence
    c = [c, dec2bin(index-1,4) - '0'];
end

end


function c = detectBPSK(r)
% Decision threshold for BPSK demodulation
threshold = 0;

% Perform detection
% If received symbol is greater than threshold, output 1; otherwise, output 0
c = (1:length(r));
for i = 1:length(r)
    if (r(i) > threshold)
        c(i) = 1;
    else
        c(i) = 0;
    end
end
end

function TEST()
% Example zero transmitted signal
transmittedSignal = zeros(1,1000);

% Example noise power
noisePower = 1;

noisySignal = channel2(transmittedSignal, noisePower);

% Calculate the power of the noisy signal
noisyPower = mean(abs(noisySignal).^2);
disp(['AWGN Noisy signal power: ', num2str(noisyPower)]);
end

function TEST_16QAM(c)
Es = 1;
s = maponto16QAM(c);
r = channel2(s,0.01);
b = detect16QAM(r);

% Define 16-QAM constellation points
P = 4*Es/(sqrt(2)+2*sqrt(10)+sqrt(18));
        
Symbols = [-3+3i,-3+1i,-3-3i,-3-1i, ...
           -1+3i,-1+1i,-1-3i,-1-1i, ...
            3+3i, 3+1i, 3-3i, 3-1i, ...
            1+3i, 1+1i, 1-3i, 1-1i]*P;

% Verify Average Symbol Energy
Es = 0;
for i = 1:length(Symbols)
    Es = Es + abs(Symbols(i));
end
Es = Es/numel(Symbols);
disp(['16QAM average symbol energy: ', num2str(Es)]);

%SER
DeterminedSymbols = maponto16QAM(b);
ErrorSymbol = sum(DeterminedSymbols ~= s);
BER = ErrorSymbol/length(s);

%BER
ErrorBit = sum(c ~= b);
SER = ErrorBit/length(c);

disp(['16QAM Empirical BER: ', num2str(SER), ', Empirical BER: ' num2str(BER)]);

% Plot Symbols
figure;
xlabel('Real');
ylabel('Imaginary');
title('16QAM Symbols with Decition Bounderies');
xlim([-5*P,5*P]);
ylim([-5*P,5*P]);
xline(0, 'r--', 'LineWidth', 1);
xline(2*P, 'r--', 'LineWidth', 1);
xline(-2*P, 'r--', 'LineWidth', 1);
yline(0, 'r--', 'LineWidth', 1);
yline(2*P, 'r--', 'LineWidth', 1);
yline(-2*P, 'r--', 'LineWidth', 1);
axis square;
grid on;
hold on;
box on;
plot(real(Symbols), imag(Symbols), 'bx','LineWidth', 2);

% Plot 
figure;
hold on;
plot(real(r), imag(r), 'gx','LineWidth', 2);
plot(real(Symbols), imag(Symbols), 'bx','LineWidth', 2);
xlabel('Real');
ylabel('Imaginary');
title('Reccived Symbols with Noise Variance: 0.01 (16QAM)');
xlim([-5*P,5*P]);
ylim([-5*P,5*P]);
xline(0, 'r--', 'LineWidth', 1);
xline(2*P, 'r--', 'LineWidth', 1);
xline(-2*P, 'r--', 'LineWidth', 1);
yline(0, 'r--', 'LineWidth', 1);
yline(2*P, 'r--', 'LineWidth', 1);
yline(-2*P, 'r--', 'LineWidth', 1);
axis square;
box on;
grid on;


end

function BERCURVES(c)
Es = 1;
No = logspace(-2,0,20);
EsNo = Es./No;
EsNo_dB = 10*log10(EsNo);

%testing
BER16QAM = (1:length(No));
for i = 1:length(No)
    %
    s = maponto16QAM(c);
    r = channel2(s,No(i));
    b = detect16QAM(r);
    
    % Claculate and store BER for plotting
    BER = sum(c ~= b)/length(c);
    BER16QAM(i) = BER;
end 

BERBPSK = (1:length(No));
for i = 1:length(No)
    %
    s = mapontoBPSK(c);
    r = channel2(s,No(i));
    b = detectBPSK(r);
    
    % Claculate and store BER for plotting
    BER = sum(c ~= b)/length(c);
    BERBPSK(i) = BER;
end

BER16QAMTHEO = 3/2*erfc(sqrt(0.1*(10.^(EsNo_dB/10))));
BERBPSKTHEO = 0.5*erfc(sqrt(EsNo)); 

figure;
semilogy(EsNo_dB, BER16QAM, '.-');
hold on;
semilogy(EsNo_dB, BER16QAMTHEO, '.-');
semilogy(EsNo_dB, BERBPSK, '.-');
semilogy(EsNo_dB, BERBPSKTHEO, '.-');
title('BER vs SNR for 16QAM and BPSK');
xlabel('Signal Noise Ratio (SNR [dB])');
ylabel('Bit Error Rate (BER)');
ylim([10^-5,1]);
grid on;
legend('16QAM', 'Theoretical 16QAM', 'BPSK', 'Theoretical BPSK', 'Location', 'southwest');


end