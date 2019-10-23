
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/inverse_3_deconvolution_sparsity')

setting = 1;
switch setting
    case 1
        % difficult
        s = 3;
        sigma = .02;
    case 2
        % easy
        s = 1.2;
        sigma = .02;
end

n = 128*2;
name = 'lena';
name = 'boat';
name = 'mri';
f0 = load_image(name);
f0 = rescale(crop(f0,n));

clf;
imageplot(f0);

x = [0:n/2-1, -n/2:-1];
[Y,X] = meshgrid(x,x);
h = exp( (-X.^2-Y.^2)/(2*s^2) );
h = h/sum(h(:));

hF = real(fft2(h));

clf;
imageplot(fftshift(h), 'Filter', 1,2,1);
imageplot(fftshift(hF), 'Fourier transform', 1,2,2);

Phi = @(x)real(ifft2(fft2(x).*hF));

y0 = Phi(f0);

clf;
imageplot(f0, 'Image f0', 1,2,1);
imageplot(y0, 'Observation without noise', 1,2,2);

y = y0 + randn(n,n)*sigma;

clf;
imageplot(y0, 'Observation without noise', 1,2,1);
imageplot(clamp(y), 'Observation with noise', 1,2,2);

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

lambda = .02;

tau = 1.5;

niter = 100;

fSpars = y;

fSpars = fSpars + tau * Phi( y-Phi(fSpars) );

fSpars = SoftThreshPsi( fSpars, lambda*tau );

exo1()

%% Insert your code here.

clf;
imageplot(clamp(fSpars), ['Sparsity deconvolution, SNR=' num2str(snr(f0,fSpars),3) 'dB']);

exo2()

%% Insert your code here.

clf;
imageplot(clamp(fBestOrtho), ['Sparsity deconvolution, SNR=' num2str(snr(f0,fBestOrtho),3) 'dB']);

J = Jmax-Jmin+1;
u = [4^(-J) 4.^(-floor(J+2/3:-1/3:1)) ];

lambda = .01;

options.ti = 1; % use translation invariance
Psi = @(a)perform_wavelet_transf(a, Jmin, -1,options);
PsiS = @(f)perform_wavelet_transf(f, Jmin, +1,options);

tau = 1.5;

a = PsiS(y);

a = a + tau * PsiS( Phi( y-Phi(Psi(a)) ) );

a = SoftThresh( a, lambda*tau );

U = repmat( reshape(u,[1 1 length(u)]), [n n 1] );
Ja = sum(sum(sum( abs(a.*U) )));

exo3()

%% Insert your code here.

fTI = Psi(a);

clf;
imageplot(fTI);

exo4()

%% Insert your code here.

clf;
imageplot(clamp(fBestTI), ['Sparsity deconvolution TI, SNR=' num2str(snr(f0,fBestTI),3) 'dB']);

exo5()

%% Insert your code here.

clf;
imageplot(clamp(fBestTV), ['TV deconvolution, SNR=' num2str(snr(f0,fBestTV),3) 'dB']);
