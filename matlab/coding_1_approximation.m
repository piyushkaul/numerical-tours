
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/coding_1_approximation')

n = 512;
f = rescale( load_image('lena', n) );

clf;
imageplot(f);

fF = fft2(f)/n;

clf;
imageplot(log(1e-5+abs(fftshift(fF))));

T = .3;
c = fF .* (abs(fF)>T);

fM = real(ifft2(c)*n);

imageplot(clamp(fM));

exo1()

%% Insert your code here.

exo2()

%% Insert your code here.

exo3()

%% Insert your code here.

Jmin = 1;
options.h = compute_wavelet_filter('Daubechies',10);
fW = perform_wavortho_transf(f,Jmin,+1, options);

clf;
plot_wavelet(fW,Jmin);
title('Wavelet coefficients');

exo4()

%% Insert your code here.

exo5()

%% Insert your code here.

fC = dct2(f);

clf;
imageplot(log(1e-5+abs(fC)));

exo6()

%% Insert your code here.

exo7()

%% Insert your code here.

w = 16;

fL = zeros(n,n);

i = 5;
j = 7;

seli = (i-1)*w+1:i*w;
selj = (j-1)*w+1:j*w;
P = f(seli,selj);

fL(seli,selj) = dct2(P);

clf;
imageplot(P,'Patch',1,2,1);
imageplot(dct2(P-mean(P(:))),'DCT',1,2,2);

exo8()

%% Insert your code here.

clf;
imageplot(min(abs(fL),.005*w*w));

exo9()

%% Insert your code here.

exo10()

%% Insert your code here.

exo11()

%% Insert your code here.

n = 512;
fList(:,:,1) = rescale( load_image('regular3',n) );
fList(:,:,2) = rescale( load_image('phantom',n) );
fList(:,:,3) = rescale( load_image('lena',n) );
fList(:,:,4) = rescale( load_image('mandrill',n) );

clf;
for i=1:4
    imageplot(fList(:,:,i),'', 2,2,i);
end

exo12()

%% Insert your code here.
