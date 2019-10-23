
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/inverse_9_tomography_sobsparse')

q = 32;

n = 256;

Theta = linspace(0,pi,q+1); Theta(end) = [];
xi = zeros(n,n);
for theta = Theta
    t = linspace(-1,1,3*n)*n;
    x = round(t.*cos(theta)) + n/2+1; y = round(t.*sin(theta)) + n/2+1;
    I = find(x>0 & x<=n & y>0 & y<=n); x = x(I); y = y(I);
    xi(x+(y-1)*n) = 1;
end

xi = fftshift(xi);

clf;
imageplot(fftshift(xi));

Phi  = @(f)fft2(f).*xi / n;

name = 'mri';
f0 = load_image(name,n);
f0 = rescale(f0);

sigma = .2;

y = Phi(f0) + sigma*randn(n,n);

PhiS = @(y)real(ifft2(y.*xi))*n;

exo1()

%% Insert your code here.

x = [0:n/2-1, -n/2:-1];
[Y,X] = meshgrid(x,x);
S = (X.^2 + Y.^2)*(2/n)^2;

lambda = 2;

fSob = real( ifft2( y .* xi ./ ( xi + lambda*S) ) )*n;

clf;
imageplot(clamp(fSob), ['Sobolev inversion, SNR=' num2str(snr(f0,fSob),3) 'dB'] );

exo2()

%% Insert your code here.

clf;
imageplot(clamp(fSob), ['Sobolev inversion, SNR=' num2str(snr(f0,fSob),3) 'dB']);

SoftThresh = @(x,T)x.*max( 0, 1-T./max(abs(x),1e-10) );

clf;
T = linspace(-1,1,1000);
plot( T, SoftThresh(T,.5) );
axis('equal');

Jmax = log2(n)-1;
Jmin = Jmax-3;

options.ti = 0; % use orthogonality.
Psi = @(a)perform_wavelet_transf(a, Jmin, -1,options);
PsiS = @(f)perform_wavelet_transf(f, Jmin, +1,options);

SoftThreshPsi = @(f,T)Psi(SoftThresh(PsiS(f),T));

clf;
imageplot( clamp(SoftThreshPsi(f0,.1)) );

lambda = 0.1;

tau = 1.5;

fSpars = PhiS(y);

fSpars = fSpars + tau * PhiS( y-Phi(fSpars) );

fSpars = SoftThreshPsi( fSpars, lambda*tau );

exo3()

%% Insert your code here.

clf;
imageplot(clamp(fSpars), ['Sparsity inversion, SNR=' num2str(snr(f0,fSpars),3) 'dB']);

exo4()

%% Insert your code here.

clf;
imageplot(clamp(fBestOrtho), ['Sparsity in Orthogonal Wavelets, SNR=' num2str(snr(f0,fBestOrtho),3) 'dB']);

options.ti = 1;
Psi = @(a)perform_wavelet_transf(a, Jmin, -1,options);
PsiS = @(f)perform_wavelet_transf(f, Jmin, +1,options);

SoftThreshPsi = @(f,T)Psi(SoftThresh(PsiS(f),T));

exo5()

%% Insert your code here.

clf;
imageplot(clamp(fBestTI), ['Sparsity in TI Wavelets, SNR=' num2str(snr(f0,fBestTI),3) 'dB']);

alpha = 3;
q = 30;
t = linspace(0,1.5,n*q*10);
x = round(.5*n*t.^alpha.*cos(2*pi*q*t)) + n/2+1; 
y = round(.5*n*t.^alpha.*sin(2*pi*q*t)) + n/2+1;
I = find(x>0 & x<=n & y>0 & y<=n); x = x(I); y = y(I);
xi = zeros(n,n); xi(x+(y-1)*n) = 1;

clf;
imageplot(xi);

exo6()

%% Insert your code here.
