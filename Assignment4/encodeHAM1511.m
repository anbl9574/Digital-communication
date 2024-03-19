function c = encodeHAM1511(b)
% b måste vara en multipel av 11 bitar
M = 4; % 4 parity bits
[H,G] = hammgen(M); 

k = reshape(b,[11,length(b)/11])' % varje rad motsvarar 11 bitar av b som ska hamming kodas
p = k*G;
result = mod(p,2); % varje rad motsvarar ett kodat meddelande

c = reshape(result',1,[])';

end

