load train;
sound(y,Fs)



%%


%Testing


Vp = 4;
N = 2;
testvector = [-4 ,1 ,5]; % Your actual test vector goes here;

% Call the function with specific parameter values
[quantizedSignal,varLin,varSat,SNqR] = MyQuantizer(testvector,Vp,N);
%%
N=2;

[bitStream] = MyGraycode(quantizedSignal,Vp,N)

