function [s] = mapontoBPSK(c)
%BPSK encodes 1 bit per symbol
Es = 1;
M = 2;
A = sqrt(Es);
dmin = 2*A; %2*Asin(pi/M)
signalConstellation = [-dmin/2,dmin/2];
k = c+1;

for i = 1:length(k)
s(i,1) = signalConstellation(k(i));
end 
end

