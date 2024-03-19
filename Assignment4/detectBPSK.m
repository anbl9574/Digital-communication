function [c] = detectBPSK(r)
Es = 1;
M = 2;
A = sqrt(Es);
dmin = 2*A; %2*Asin(pi/M)
signalConstellation = [-dmin/2,dmin/2];

dist = cell(length(r),1);
dist(:) = {signalConstellation};

%eftersom signalConstellation bara har en rad fås m = 1 för alla u i detta
%fall
for u = 1:length(r)
    dist{u} = dist{u}-r(u);
    minimum(u)=min(min(abs(dist{u}))); % find smallest distance to symbol in constellation
    [m(u,1),n(u,1)]= find(abs(dist{u}) == minimum(u)); %find index of smallest distance to symbol in constellation
end

bin = [0,1];
for i = 1:length(n)
c(i,1) = bin(n(i));
end 
end

