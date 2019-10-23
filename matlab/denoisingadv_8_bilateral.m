
warning off
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/denoisingadv_8_bilateral')
warning on

n = 256*2;

name = 'hibiscus';
f0 = load_image(name, n);
f0 = rescale(crop( sum(f0,3) ,n));

x = [0:n/2-1, -n/2:-1];
[Y,X] = meshgrid(x,x);
GaussianFilt = @(s)exp( (-X.^2-Y.^2)/(2*s^2) );

clf;
imageplot(fftshift(GaussianFilt(40)));

Filter = @(F,s)real( ifft2( fft2(F).*repmat( fft2(GaussianFilt(s)), [1 1 size(F,3)] ) ) );

clf;
imageplot( Filter(f0,5) );

exo1()

%% Insert your code here.

sx = 5;

sv = .2;

p = 10;

Gaussian = @(x,sigma)exp( -x.^2 / (2*sigma^2) );
WeightMap = @(f0,sv)Gaussian( repmat(f0, [1 1 p]) - repmat( reshape((0:p-1)/(p-1), [1 1 p]) , [n n 1]), sv );

W = WeightMap(f0,sv);

exo2()

%% Insert your code here.

bileteral_stack_tmp = @(f0,sx,W)Filter(W.*repmat(f0, [1 1 p]), sx) ./ Filter(W, sx);
bileteral_stack = @(f0,sx,sv)bileteral_stack_tmp(f0,sx,WeightMap(f0,sv));

F = bileteral_stack(f0,sx,sv);

exo3()

%% Insert your code here.

[y,x] = meshgrid(1:n,1:n);
indexing = @(F,I)F(I);
destacking = @(F,I)indexing(F,x + (y-1)*n + (I-1)*n^2);

bilateral_nn = @(f0,sx,sv)destacking( bileteral_stack(f0,sx,sv), round( f0*(p-1) ) + 1 );

fNN = bilateral_nn(f0,sx,sv);
clf;
imageplot( fNN );

frac = @(x)x-floor(x);
lininterp = @(f1,f2,Fr)f1.*(1-Fr) + f2.*Fr;
bilateral_lin1 = @(F,f0)lininterp( destacking(F, clamp(floor(f0*(p-1)) + 1,1,p) ), ...
                                  destacking(F, clamp(ceil(f0*(p-1)) + 1,1,p) ), ...
                                  frac(f0*(p-1)) );
bilateral_lin = @(f0,sx,sv)bilateral_lin1(bileteral_stack(f0,sx,sv), f0);

exo4()

%% Insert your code here.

exo5()

%% Insert your code here.

exo6()

%% Insert your code here.

mu = .05;

f = f0 + randn(n,n)*mu;

clf;
imageplot(clamp(f));

sx = 4; sv = .2;
clf;
imageplot( clamp(bilateral_lin(f,sx,sv)) );

exo7()

%% Insert your code here.

clf;
imageplot(clamp(fOpt), ['Bilateral, SNR=' num2str(snr(f0,fOpt),3) 'dB']);

exo8()

%% Insert your code here.

sx = 4;
sv = .2;

f1 = bilateral_lin(f0,sx,sv);

r = f0 - f1;

clf;
imageplot(f1);

clf;
imageplot(r);

exo9()

%% Insert your code here.

exo10()

%% Insert your code here.

addpath('toolbox_additional/');
name = 'memorial';
f = load_hdr([name '.hdr']);
p = min(size(f,1),size(f,2));
f = rescale( f(1:p,1:p,:) );
f = rescale( clamp( image_resize(f,n,n) ) );

clf;
imageplot(min(f,1e-4));

clf;
imageplot(min(f,1e-3));

fhsv = rgb2hsv(f);
fV = fhsv(:,:,3);

color_recompose = @(fV)hsv2rgb( cat(3, fhsv(:,:,1:2), rescale(fV) ) );

epsilon = 1e-5;
alpha = 0.1;
imageplot( color_recompose( (fV+epsilon).^alpha ) );

exo11()

%% Insert your code here.

epsilon = 1e-5;
FV = log(fV+epsilon);
a = min(FV(:)); b = max(FV(:));
FV = (FV-a)/(b-a);

sx = 5; sv = .02;
FV1 = bilateral_lin(FV, sx,sv);

clf;
imageplot(FV1);

clf;
imageplot(FV-FV1);

exo12()

%% Insert your code here.
