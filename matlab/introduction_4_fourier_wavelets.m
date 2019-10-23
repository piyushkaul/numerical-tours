
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/introduction_4_fourier_wavelets')

name = 'lena';
n0 = 512;
f = rescale( load_image(name,n0) );

clf;
imageplot( f, 'Image f');

clf;
imageplot( crop(f,64), 'Zoom' );

clf;
imageplot(-f, '-f', 1,2,1);
imageplot(f(n0:-1:1,:), 'Flipped', 1,2,2);

k = 9; % size of the kernel
h = ones(k,k);
h = h/sum(h(:)); % normalize

fh = perform_convolution(f,h);

clf;
imageplot(fh, 'Blurred image');

F = fft2(f) / n0;

disp(strcat(['Energy of Image:   ' num2str(norm(f(:)))]));
disp(strcat(['Energy of Fourier: ' num2str(norm(F(:)))]));

L = fftshift(log( abs(F)+1e-1 ));

clf;
imageplot(L, 'Log(Fourier transform)');

M = n0^2/64;

exo1()

%% Insert your code here.

clf;
subplot(2,1,1);
plot(f(:,n0/2)); 
axis('tight'); title('f');
subplot(2,1,2);
plot(fM(:,n0/2)); 
axis('tight'); title('f_M');

T = .2;

F = fft2(f) / n0;

FT = F .* (abs(F)>T);

clf;
L = fftshift(log( abs(FT)+1e-1 ));
imageplot(L, 'thresholded Log(Fourier transform)');

fM = real( ifft2(FT)*n0 );

clf;
imageplot(clamp(fM), ['Non-linear, Fourier, SNR=' num2str(snr(f,fM), 4) 'dB']);

m = sum(FT(:)~=0);
disp(['M/N = 1/'  num2str(round(n0^2/m)) '.']);

exo2()

%% Insert your code here.

Jmin = 0;

fw = perform_wavelet_transf(f,Jmin,+1);

clf;
plot_wavelet(fw);

exo3()

%% Insert your code here.

T = .2;

fwT = fw .* (abs(fw)>T);

clf;
subplot(1,2,1);
plot_wavelet(fw);
title('Original coefficients');
subplot(1,2,2);
plot_wavelet(fwT);

fM = perform_wavelet_transf(fwT,Jmin,-1);

clf;
imageplot(clamp(fM), strcat(['Approximation, SNR=' num2str(snr(f,fM),3) 'dB']));

exo4()

%% Insert your code here.

clf;
subplot(2,1,1);
plot(f(:,n0/2)); 
axis('tight'); title('f');
subplot(2,1,2);
plot(fM(:,n0/2)); 
axis('tight'); title('f_M');
