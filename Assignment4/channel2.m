function [noisySignal] = channel2(transmittedSignal,noisePower)
%observera att variansen av AGWN är lika med dess effekt.
N = length(transmittedSignal);
Pw = sqrt(noisePower).*(randn(1,N)+ j*randn(1,N))./sqrt(2);%se fotnot 3 i A3
if iscolumn(transmittedSignal) == 1
    Pw  = transpose(Pw);
end
noisySignal = transmittedSignal + Pw;
end
