clear all
close all
clc
%% Reading the EEG signal data from the csv file with 14 columns of all data channels

EEG_data_C3_4 = csvread('test.csv');
V = EEG_data_C3_4(:,1);
N=length(EEG_data_C3_4(:,1));

%% Normalising the EEG Data Signal
% Normalized EEG signal = Collected EEG Signal – Mean (Collected EEGSignal)
%                               Std (Collected EEG Signal)
m1 = mean(V);
std1 = std(V); 
Normalised_V = (V-m1)/std1;
figure;
subplot(2,1,1);plot(V);title('Original'); 
subplot(2,1,2);plot(Normalised_V);title('Normalised Signal'); 

%% Perfmorming De-noising

% Plotting frequency spectrum of our original signal
Fs =1000;
X = fft(V);
P2 = abs(X/N); % two sided spectrum
P1 = P2(1:N/2+1); % single sided spectrum
P1(2:end-1) = 2*P1(2:end-1)
f = Fs*(0:(N/2))/N;
figure;
subplot(2,1,1)
plot(f, P1)
xlabel('f (Hz)')
ylabel('|V(f)|')
title('Single-Sided Amplitude Spectrum of V(t)')

%Applying bandpass filter to filter out the unwanted signal <4 and >30Hz
yy1=bandpass(V,[4 30],128);
waveletFunction = 'db8';

% Plotting frequency spectrum of the denoised signal

X1 = fft(yy1);
P21 = abs(X1/N); % two sided spectrum
P11 = P21(1:N/2+1); % single sided spectrum
P11(2:end-1) = 2*P11(2:end-1)

f1 = Fs*(0:(N/2))/N;
subplot(2,1,2)
plot(f1, P11)
xlabel('f (Hz)')
ylabel('|V_(f)|')
title('Single-Sided Amplitude Spectrum of denoised V(t)')


%% Performing Discrete Wavelet Transform
% Perfmorming DISCRETE WAVELET TRANSFORM by using a 8 levelwavelet
[C,L] = wavedec(V,8,waveletFunction); % Denoise EOG from EEG
% we obtain the Approximation(A) & Detail(D) coefficients at various levels
cD1 = detcoef(C,L,1);
cD2 = detcoef(C,L,2);
cD3 = detcoef(C,L,3);
cD4 = detcoef(C,L,4);
cD5 = detcoef(C,L,5); %GAMMA
cD6 = detcoef(C,L,6); %BETA
cD7 = detcoef(C,L,7); %ALPHA
cD8 = detcoef(C,L,8); %THETA
cA8 = appcoef(C,L,waveletFunction,5); %DELTA

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
figure;
subplot(5,2,[1 2]);plot(V);title('Original'); 
subplot(5,2,3);plot(A1);title('Level 1 Approximation');
subplot(5,2,4);plot(A2);title('Level 2 Approximation');
subplot(5,2,5);plot(A3);title('Level 3 Approximation');
subplot(5,2,6);plot(A4);title('Level 4 Approximation');
subplot(5,2,7);plot(A5);title('Level 5 Approximation');
subplot(5,2,8);plot(A6);title('Level 6 Approximation');
subplot(5,2,9);plot(A7);title('Level 7 Approximation');
subplot(5,2,10);plot(A8);title('Level 8 Approximation');
%axis off
sgtitle('wrcoef Approximation vs Original eeg');

%% Plotting different bands of frequency of input signals:

figure;
subplot(6,1,1);plot(V);title('Original');
subplot(6,1,2);plot(D5);title('Gamma');
subplot(6,1,3);plot(D6);title('Beta');
subplot(6,1,4);plot(D7);title('Alpha');
subplot(6,1,5);plot(D8);title('Theta');
subplot(6,1,6);plot(A8);title('Delta');
BetaSize = size(D6)

%%
fs=128;
fmin = 30; % minimum passband frequency in Hz (High Gamma) 30
fmax =50; % maximum passband frequency in Hz (High Gamma) 50
Rs = 20; % stopband attenuation in dB 20
Rp = 1; % passband ripples in dB 1
% for High gamma band
[order, Wn] = ellipord([fmin/(fs/2),fmax/(fs/2)],[(fmin-1)/(fs/2),(fmax+1)/(fs/2)],Rp,Rs);
[B,A] = ellip(order, Rp, Rs, [fmin/(fs/2), fmax/(fs/2)]);
%filtering of entire data into High Gamma band
data_hgamma = filtfilt(B,A,V);
figure;
subplot(4,1,1);plot(D6);title('wavelet');%axis off
subplot(4,1,2);plot(data_hgamma);title('ellipord');%axis off
sgtitle('waverec gamma vs ellipord gamma');
dLen = length(D7);
dF = fs/dLen;
f = dF*(0:dLen-1)';
A0 = waverec(C,L,'db8');
err = max(abs(V-A0));
disp(err);
Cnew = C;
tmp = cumsum(L);
Cnew(tmp(end-2)+1:tmp(end-1)) = 0;
Rec_signal=waverec(Cnew,L,'db8');
subplot(4,1,3);plot(V,'k');title('Original'); %axis off
subplot(4,1,4); plot(Rec_signal,'r');title('Reconstruct');%axis off
%subplot(3,1,3);plot(Rec_signal);title('Rec_signal');%axis off
% figure
% subplot(4,1,2); plot(Rec_signal,'r','linewidth',2);title('Reconstruct');axis off
% sgtitle('waverec reconstruct wavelet');



