function c = encodeHAM74(b)
% b and c are column vectors!
% b must be evenly divisible by 2!
%
%Theory----------
%https://www.ece.unb.ca/tervo/ee4253/hamming2.shtml
%https://www.ece.unb.ca/tervo/ee4253/hamming3.shtml
% n = 2^M-1, and k = n - M. For example, M=3 produces a (7,4) code.
%
%Misc------------
%https://www.csus.edu/indiv/p/pangj/166/f/d8/5_Modulo%202%20Arithmetic.pdf
%https://en.wikipedia.org/wiki/Finite_field_arithmetic

M = 3; % 3 parity bits
[H,G] = hammgen(M); 

k = reshape(b,[4,length(b)/4])' % varje rad motsvarar 4 bitar av b som ska hamming kodas
p = k*G;
result = mod(p,2); % varje rad motsvarar ett kodat meddelande

c = reshape(result',1,[])';

%jämför med CODE = encode(b, 7, 4, 'hamming')
%sum(c-CODE) ska ge 0 om alla lika

end

