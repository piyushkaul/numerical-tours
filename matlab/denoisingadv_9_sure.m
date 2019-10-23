
warning('off')
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/denoisingadv_9_sure')
warning('on')

n = 128*2;
N = n^2;

name = 'hibiscus';
f0 = rescale( sum(load_image(name,n),3) );

clf;
imageplot(f0);

sigma = .08;

f = f0 + sigma*randn(n);

clf;
imageplot(clamp(f), strcat(['Noisy, SNR=' num2str(snr(f0,f),3) 'dB']));

convol = @(f,g)real(ifft2(fft2(f) .* repmat(fft2(g), [1 1 size(f,3)]) ));

normalize = @(f)f/sum(f(:));
x = [0:n/2 -n/2+1:-1];
[Y,X] = meshgrid(x,x);
g = @(lambda)normalize( exp( -(X.^2+Y.^2)/(2*lambda^2) ) );

h = @(f,lambda)convol(f, g(lambda));

lambda = 1.5;
clf;
imageplot(clamp(h(f,lambda)));

df = @(lambda)real(sum(sum(fft2(g(lambda)))));

SURE = @(f,hf,lambda)-N*sigma^2 + norm(hf-f, 'fro')^2 + 2 * sigma^2 * df(lambda);

exo1()

%% Insert your code here.

exo2()

%% Insert your code here.

exo3()

%% Insert your code here.

f = f0 + sigma*randn(n);

W  = @(f)perform_wavortho_transf(f,0,+1);
Ws = @(x)perform_wavortho_transf(x,0,-1);

clf;
plot_wavelet(W(f0),1);

S = @(x,lambda)max(0, 1-lambda ./ max(1e-9,abs(x)) ) .* x;

h = @(f,lambda)Ws(S(W(f),lambda));

lambda = 3*sigma/2;
clf;
imageplot(clamp(h(f,lambda)));

df = @(hf,lambda)sum(sum( abs(W(hf))>1e-8 ));

SURE = @(f,hf,lambda)-N*sigma^2 + norm(hf-f, 'fro')^2 + 2 * sigma^2 * df(hf,lambda);

exo4()

%% Insert your code here.

exo5()

%% Insert your code here.

exo6()

%% Insert your code here.

q = 4;

[dX,dY,X,Y] = ndgrid(0:q-1,0:q-1,1:q:n-q+1,1:q:n-q+1);
I = X+dX + (Y+dY-1)*n;
blocks = @(fw)reshape(fw(I(:)),size(I));

linearize = @(x)x(:);
unblock = @(H)reshape( accumarray( I(:), linearize(H), [n*n 1], @min), [n n]);

energy = @(H)mean(mean(abs(H).^2,1),2);
energy = @(H)repmat( max3(energy(H),1e-15), [q q]);

S = @(H,lambda)max(1-lambda^2 ./ energy(H),0) .* H;

h = @(f,lambda)Ws(unblock(S(blocks(W(f)),lambda) ) );

lambda = 1.1*sigma;
clf;
imageplot(clamp(h(f,lambda)));


