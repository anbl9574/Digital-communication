function bEst = decodeHAM74(c)
% b and c are column vectors!
% c must be evenly divisible by 2!

M = 3; % 3 parity bits
[H,G] = hammgen(M); 

%OBS!
% Olika varianter på Hamming-kod finns (olika generator matriser) https://en.wikipedia.org/wiki/Hamming(7,4) vilket gör att 
% datan och paritetsbitarna hamnar på olika ställen i de kodade orden:
% wiki: p1 p2 d1 p3 d2 d3 d4
% matlab hammgen: p1 p2 p3 d1 d2 d3 d4 (bra att veta när man avkodar...)

%letar upp och korrigerar fel innan avkodning (kan max korrigera 1 bitfel i varje meddelande): 

rM = reshape(c,[7,length(c)/7]); %varje column motsvarar ett kodat meddelande 7 rader * # kolumner

s1 = H*rM; %beräknar syndromet (utan modulo 2)

s = mod(s1,2); % varje kolumn i syndromet som inte är noll motsvarar ett fel i det korresponderande meddelandet
               % Kolumnens nummer i H motsvarar platsen felet finns på i rM
               % (för varje meddelande alltså) 
           
[~,ind] = ismember(s',H','rows');

for j = 1:length(ind)  
  if ind(j) ~= 0
    rM(ind(j),j) = not(rM(ind(j),j)); %fixar felet.
  end
end

bEst = [];
for l = 1:size(rM,2) %för varje rad i rM
    bEst = vertcat(bEst,rM(4:end,l)); %plockar ut de sista 4a elementen.
end 
%bEst = bEst';

end

