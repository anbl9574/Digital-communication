function[bitStream] = MyGraycode(quantizedSignal,Vp,N)

step=(Vp*2)/2^N;
maxlvl = (2^N/2-0.5)*step;
minlvl = -(2^N/2-0.5)*step;

levelquant = linspace(minlvl,maxlvl,2^N);

for i=1:length(levelquant)
    quantizedSignal(    


num = dec2bin(quantizedSignal, N) - '0'; % Convert to numeric array

bitStream = bitxor(num, bitshift(num, -1));

