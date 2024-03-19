
%enkelt testfall
% x = [1 0 1 1 0 0 1 0 1 1 1 1]
% 
% x = randi([0 1],100,1); %jämnt delbart på 4 ---> använd för Ham74
% x = randi([0 1],55,1); %jämnt delbart på 11 ---> använd för Ham1511


%BER vs SNR curves
%Från uppgiftsbeskrivningen ska två teoretiska kurvor plottas för 
%BPSK + HAM74 och BPSK + HAM1511
%Eb_No ska innehålla jämnt fördelade (logaritmiska) datapunkter mellan
% A = 10log10(P2/P1)
% a/b = log10(dB)
%-10dB (0.316228) och 15dB (31.6228) ==> a =-0.5  b = 1.5 i logspace()

Eb_No = logspace(-0.5,1.5,1000);
%Eb_No_dB = 10*log10(Eb_No);
n = 7;
k = 4;
for i = 1:length(Eb_No)
    BER_BPSK_HAM74(1,i) = 0.5*erfc(sqrt((k/n)*Eb_No(i)))*(1 - (1-0.5*erfc(sqrt((k/n)*Eb_No(i))))^(n-1));
end

n = 15;
k = 11;
for i = 1:length(Eb_No)
    BER_BPSK_HAM1511(1,i) =  0.5*erfc(sqrt((k/n)*Eb_No(i)))*(1 - (1-0.5*erfc(sqrt((k/n)*Eb_No(i))))^(n-1));
end



semilogy(db(Eb_No,'power'),BER_BPSK_HAM74);
ylim([10^-7,10^-0.7014])
title('Theoretical BER for BPSK using 74-HAM and 1511-HAM')
ylabel('BER')
xlabel('Eb/No (db)')
hold on
semilogy(db(Eb_No,'power'),BER_BPSK_HAM1511);
legend('74-HAM','1511-HAM')




%bitström
%encodeham
%maptobpsk
%channel(var)
%detectBpsk
%decodeham


%skräp

% k = reshape(x,[4,length(x)/4])
% k = reshape(x,[],)
% A = [1 2 ; 3 4]
% reshape(A,1,[])
% transpose(A)
% A.'
% A(:)
% reshape(A.',1,[])


%%

n = 7;
k = 4;
% length is a multiple of 4 ( HAM74 )
x = randi ([0 1] ,1000 ,1) ;
c = encodeHAM74 ( x ) ;
s = mapontoBPSK ( c ) ;

Eb_No = logspace ( -0.5 ,1.5 ,1000) ;
Eb = 1;
No = Eb ./ Eb_No ;
noiseVar = No ;

for l = 1: length ( noiseVar )
    nosiySignal (: , l ) = channel2 (s , noiseVar ( l ) ) ;
    c_detect (: , l ) = detectBPSK ( nosiySignal (: , l ) ) ;
    bEst (: , l ) = decodeHAM74 ( c_detect (: , l ) ) ;
    bitErrors (l ,1) = sum ( abs (x - bEst (: , l ) ) ) ;
    BER (l ,1) = bitErrors (l ,1) / length ( x ) ;
end

for i = 1: length ( Eb_No )
    BER_BPSK_HAM74 (1 , i ) = 0.5* erfc ( sqrt (( k / n ) * Eb_No ( i )) ) *(1 - (1 -0.5* erfc ( sqrt (( k / n ) * Eb_No ( i ) ) ) ) ^( n-1) ) ;
end


semilogy(db(Eb_No ,'power'),BER_BPSK_HAM1511);
hold on
semilogy(db(Eb_No ,'power'),BER_BPSK_HAM74);
ylim ([10^ -7 ,10^ -0.7014])
title (' Theoretical and Empirical BER for BPSK using - HAM ')
ylabel ('BER ')
xlabel ('Eb/No (db)')
hold on
semilogy(db(Eb_No, 'power'), BER, '.');
legend (' Theoretical BER ','Empirical BER ') ;