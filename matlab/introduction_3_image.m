
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/introduction_3_image')

name = 'lena';
n = 256;
M = load_image(name, []);
M = rescale(crop(M,n));

clf;
imageplot(M, 'Original', 1,2,1);
imageplot(crop(M,50), 'Zoom', 1,2,2);

clf;
imageplot(-M, '-M', 1,2,1);
imageplot(M(n:-1:1,:), 'Flipped', 1,2,2);

k = 9;
h = ones(k,k);
h = h/sum(h(:));

Mh = perform_convolution(M,h);

clf;
imageplot(M, 'Image', 1,2,1);
imageplot(Mh, 'Blurred', 1,2,2);

G = grad(M);
clf;
imageplot(G(:,:,1), 'd/dx', 1,2,1);
imageplot(G(:,:,2), 'd/dy', 1,2,2);

Mf = fft2(M);
Lf = fftshift(log( abs(Mf)+1e-1 ));
clf;
imageplot(M, 'Image', 1,2,1);
imageplot(Lf, 'Fourier transform', 1,2,2);

exo1()

%% Insert your code here.

exo2()

%% Insert your code here.

p = 64;
n = p*4;
M = load_image('boat', 2*p); M = crop(M,p);
Mf = fftshift(fft2(M));
MF = zeros(n,n);
sel = n/2-p/2+1:n/2+p/2;
sel = sel;
MF(sel, sel) = Mf;
MF = fftshift(MF);
Mpad = real(ifft2(MF));
clf;
imageplot( crop(M), 'Image', 1,2,1);
imageplot( crop(Mpad), 'Interpolated', 1,2,2);

Mspline = image_resize(M,n,n);
clf;
imageplot( crop(Mpad), 'Fourier (sinc)', 1,2,1);
imageplot( crop(Mspline), 'Spline', 1,2,2);
