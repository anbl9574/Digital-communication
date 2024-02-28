close all
clearvars
clc

transmittedSignal = zeros(1, 1000); % Example zero transmitted signal
symbolPower = 1; % Example noise power

noisySignal = channel2(transmittedSignal, symbolPower);

% Calculate the power of the noisy signal
noisyPower = mean(abs(noisySignal).^2);

disp(['Noisy signal power: ', num2str(noisyPower)]);

%QPSK

% Number of symbols
numSymbols = 4;
symbol = 01 + 1;

% Specify phase angles for QPSK symbols
phase_angles = [pi/4, 3*pi/4, 5*pi/4, 7*pi/4]; % Phase angles for QPSK symbols

% Randomly choose symbols
for i = 1:numSymbols
    symbol_indices(i) = i;
end

% Create complex symbols with unit magnitude and specified phases
symbols = exp(1i * phase_angles(symbol_indices));

% Verify average symbol power
average_power = mean(abs(symbols).^2);
disp(['Average symbol power: ', num2str(average_power)]);

% Plot the symbols and detection regions in the complex plane
figure;
plot(real(symbols), imag(symbols), 'rx','LineWidth', 2);
xlabel('In-phase');
ylabel('Quadrature');
title('QPSK Symbols with Detection Regions');
xlim([-1,1]);
ylim([-1,1]);
axis square;
grid on;
hold on;

% Define detection regions and plot
xline(0, 'r--', 'LineWidth', 2);
yline(0, 'r--', 'LineWidth', 2);
% for i = 1:length(phase_angles)/2
%     % Calculate the midpoint between adjacent symbols
%     mid_angle = (phase_angles(i) + phase_angles(mod(i, 4) + 1)) / 2;
%     
%     % Calculate the vertices of the square detection region
%     x = [-2*sqrt(symbolPower / 2), 2*sqrt(symbolPower / 2), 2*sqrt(symbolPower / 2), -2*sqrt(symbolPower / 2)] * cos(mid_angle);
%     y = [-2*sqrt(symbolPower / 2), -2*sqrt(symbolPower / 2), 2*sqrt(symbolPower / 2), 2*sqrt(symbolPower / 2)] * sin(mid_angle);
%     
%     % Plot the detection region
%     plot(x, y, 'r--','LineWidth', 2);
% end
hold off;

SignalSymbols = symbols(2) * ones(1, 1000);
noisePower2 = [0.01, 0.1, 1];
for i = 1:length(noisePower2)
    ReccivedSymbols = channel2(SignalSymbols,noisePower2(i));
    figure;
    scatter(real(ReccivedSymbols), imag(ReccivedSymbols), 'b.');
    hold on;
    plot(real(SignalSymbols), imag(SignalSymbols), 'rx', 'MarkerSize', 5, 'LineWidth', 2);
    xlabel('Real-part');
    ylabel('Imaginary-part');
    
    title('QPSK Sent Vs Reccived Symbols with Detection Regions');
    xline(0, 'r--', 'LineWidth', 2);
    yline(0, 'r--', 'LineWidth', 2);

    xmax = max(real(ReccivedSymbols));
    xmin = min(real(ReccivedSymbols));
    xlims = max(abs(xmax),abs(xmin));
    lim = xlims + 0.1;
    xlim([-lim, lim]);
    ylim([-lim, lim]);
    
    legend('Noisy Reccived Symbols', 'Sent Symbol' ,'Decition Boundries')
    axis square;
    grid on;
    hold off;
    
    %SER
    ErrorSymbol = sum(angle(ReccivedSymbols) > pi | angle(ReccivedSymbols) < pi/2);
    SER = ErrorSymbol/length(ReccivedSymbols);
    disp(['Noise varriance: (', num2str(noisePower2(i)),') Empirical SER: ',num2str(SER)]);
end


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