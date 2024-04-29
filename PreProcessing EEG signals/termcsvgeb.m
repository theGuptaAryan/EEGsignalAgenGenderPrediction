EEG_data_C3_4 = csvread('DATASET_ELL319.csv');
N=length(EEG_data_C3_4(:,i));
ma = zeros(N,14);
for i = 1:14    
    V = EEG_data_C3_4(:,i);
    N=length(EEG_data_C3_4(:,i));
    yy1=bandpass(V,[4 30],128);
    waveletFunction = 'db8';
    [C,L] = wavedec(V,8,waveletFunction);
    % Reconstruct single branch from 1-D wavelet coefficients
    D1 = wrcoef('d',C,L,waveletFunction,1); %NOISY
    A1 = wrcoef('a',C,L,waveletFunction,1); %Level 1 Approximation
    D2 = wrcoef('d',C,L,waveletFunction,2); %NOISY
    A2 = wrcoef('a',C,L,waveletFunction,2); %Level 2 Approximation
    D3 = wrcoef('d',C,L,waveletFunction,3); %NOISY
    A3 = wrcoef('a',C,L,waveletFunction,3); %Level 3 Approximation
    D4 = wrcoef('d',C,L,waveletFunction,4); %NOISY
    A4 = wrcoef('a',C,L,waveletFunction,4); %Level 4 Approximation
    D5 = wrcoef('d',C,L,waveletFunction,5); %GAMMA
    A5 = wrcoef('a',C,L,waveletFunction,5); %Level 5 Approximation
    D6 = wrcoef('d',C,L,waveletFunction,6); %BETA
    A6 = wrcoef('a',C,L,waveletFunction,6); %Level 6 Approximation
    D7 = wrcoef('d',C,L,waveletFunction,7); %ALPHA
    A7 = wrcoef('a',C,L,waveletFunction,7); %Level 7 Approximation
    D8 = wrcoef('d',C,L,waveletFunction,8); %THETA 4-7
    A8 = wrcoef('a',C,L,waveletFunction,8); %DELTA 1-4 Level 8 Approximation
    ma(:,i) = D6;
end
csvwrite('all.csv',ma)
111