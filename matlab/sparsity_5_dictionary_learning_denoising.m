
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/sparsity_5_dictionary_learning_denoising')

sigma = .06;

n0 = 256;

name = 'barb';
f0 = rescale( crop(load_image(name),n0) );

clf;
imageplot(f0);

f = f0 + sigma*randn(n0);

clf;
imageplot(f);

if not(exist('D'))
    addpath('m_files/')
    sparsity_4_dictionary_learning;
end

w = sqrt(size(D,1));

p = size(D,2);

n = w*w;

q = 2;

[y,x] = meshgrid(1:q:n0-w/2, 1:q:n0-w/2);
m = size(x(:),1);
Xp = repmat(dX,[1 1 m]) + repmat( reshape(x(:),[1 1 m]), [w w 1]);
Yp = repmat(dY,[1 1 m]) + repmat( reshape(y(:),[1 1 m]), [w w 1]);

Xp(Xp>n0) = 2*n0-Xp(Xp>n0);
Yp(Yp>n0) = 2*n0-Yp(Yp>n0);

Y = f(Xp+(Yp-1)*n0);
Y = reshape(Y, [n, m]);

theta = mean(Y);
Y = Y - repmat( theta, [n 1] );

rho = .95;
epsilon = rho*sqrt(n)*sigma;

U = (eye(p) + D'*D)^(-1);
Replicate = @(z)deal(z, D*z);
ProjC = @(x,u)Replicate( U*( x + D'*u ) );

ProxG = @(f,u,gamma)ProjC(f,u);

y = Y(:,1:m);

amplitude = @(a)repmat( sqrt(sum(a.^2,1)), [n 1] );
ProjB = @(u)y + (u-y) .* min(1, epsilon./amplitude(u-y) );

ProxL1 = @(x,gamma)max(0,1-gamma./max(1e-9, abs(x))) .* x;

ProxF = @(x,u,gamma)deal( ProxL1(x,gamma), ProjB(u) );

mu = 1;
gamma = 1;

niter = 800;

exo1()

%% Insert your code here.

Y1 = reshape(D*x, [w w m]);

Y1 = Y1 - repmat( mean(mean(Y1)), [w w] );
Y1 = Y1 + reshape(repmat( theta, [n 1] ), [w w m]);

W = zeros(n0,n0);
f1 = zeros(n0,n0);
for i=1:m
    x = Xp(:,:,i); y = Yp(:,:,i);
    f1(x+(y-1)*n0) = f1(x+(y-1)*n0) + Y1(:,:,i);
    W(x+(y-1)*n0) = W(x+(y-1)*n0) + 1;
end
f1 = f1 ./ W;

clf;
imageplot(clamp(f1), ['Denoised, SNR=' num2str(snr(f0,f1),4) 'dB']);

exo2()

%% Insert your code here.

exo3()

%% Insert your code here.
