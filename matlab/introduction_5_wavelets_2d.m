
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/introduction_5_wavelets_2d')

name = 'cortex';
n0 = 512;
f = load_image(name,n0);
f = rescale( sum(f,3) );

clf;
imageplot(f);

Jmin = 0;

Psi = @(f)perform_wavelet_transf(f,Jmin,+1);

PsiS = @(fw)perform_wavelet_transf(fw,Jmin,-1);

fW = Psi(f);

clf;
plot_wavelet(fW);

T = .5;

Thresh = @(fW,T)fW .* (abs(fW)>T);

fWT = Thresh(fW,T);

exo1()

%% Insert your code here.

clf;
subplot(1,2,1);
plot_wavelet(fW);
title('Original coefficients');
subplot(1,2,2);
plot_wavelet(fWT);

f1 = PsiS(fWT);

clf;
imageplot(f, 'Image', 1,2,1);
imageplot(clamp(f1), strcat(['Approximation, SNR=' num2str(snr(f,f1),3) 'dB']), 1,2,2);

M = n0^2/16;

exo2()

%% Insert your code here.

fWT = Thresh(fW,T);

disp(strcat(['      M=' num2str(M)]));
disp(strcat(['|fWT|_0=' num2str(sum(fWT(:)~=0))]));

exo3()

%% Insert your code here.

sigma = .1;

y = f + randn(n0,n0)*sigma;

clf;
imageplot(f, 'Clean image', 1,2,1);
imageplot(clamp(y), ['Noisy image, SNR=' num2str(snr(f,y),3) 'dB'], 1,2,2);

fW = Psi(y);

T = 3*sigma;

fWT = Thresh(fW,T);

clf;
subplot(1,2,1);
plot_wavelet(fW);
title('Original coefficients');
subplot(1,2,2);
plot_wavelet(fWT);

f1 = PsiS(fWT);

clf;
imageplot(clamp(y), 'Noisy image', 1,2,1);
imageplot(clamp(f1), strcat(['Denoising, SNR=' num2str(snr(f,f1),3) 'dB']), 1,2,2);

exo4()

%% Insert your code here.

tau_max = 8;

[Y,X] = meshgrid(0:tau_max-1,0:tau_max-1);

f1 = zeros(n0,n0);

i = 1;

fTrans = circshift(y,[X(i) Y(i)]);

fTrans = PsiS( Thresh( Psi(fTrans) ,T) );

fTrans = circshift(fTrans,-[X(i) Y(i)]);

f1 = (i-1)/i*f1 + fTrans/i;

exo5()

%% Insert your code here.

exo6()

%% Insert your code here.

exo7()

%% Insert your code here.
