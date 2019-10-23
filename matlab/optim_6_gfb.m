
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/optim_6_gfb')

name = 'lena';
N = 256;
f0 = load_image(name);
f0 = rescale(crop(f0,N));

clf
imageplot(f0);

rho_M = .7;
mask = rand(N,N) > rho_M;
M = @(f) mask.*f;

sig_K = 2;
[X,Y] = meshgrid( [0:N/2-1 -N/2:-1] );
k = exp( - (X.^2+Y.^2) / (2*sig_K^2) );
k = k./sum(abs(k(:)));

Fk = fft2(k);
K = @(f) real( ifft2( fft2(f).*Fk ) );

Phi = @(f)M(K(f));

sig_w = .025;
y = Phi(f0) + sig_w*randn(N,N);

clf
imageplot(y);

G = @(u) sqrt( sum( u.^2, 3 ) );
proxG = @(u,gamma) repmat( max(0,1 - gamma./G(u) ), [1 1 2] ).*u;

L = @(x) cat( 3, x(2:2:end,2:2:end)-x(1:2:end,1:2:end), x(1:2:end,2:2:end)-x(2:2:end,1:2:end) );

LShift = @(x,s) L( circshift(x,s) );
Li = {@(x)LShift(x,[0,0]), @(x)LShift(x,[1,0]), @(x)LShift(x,[0,1]), @(x)LShift(x,[1,1]) };

for i=1:4
    Gi{i} = @(x) sum( sum( G( Li{i}(x) ) ) );
end

U = @(x)upsampling( upsampling( x, 1, 2 ), 2, 2 );

revIdx = [N 1:N-1];
gradS = @(u)(u(revIdx,revIdx,1) - u(:,:,1)) + (u(:,revIdx,2) - u(revIdx,:,2));

L1S = @(v)gradS(U(v));

LShiftS = @(gx,s) circshift( L1S( gx ), -s ); 
LiS = { @(x)LShiftS(x,[0,0]), @(x)LShiftS(x,[1,0]), @(x)LShiftS(x,[0,1]), @(x)LShiftS(x,[1,1]) };

exo1()

%% Insert your code here.

proxG_Id = @(u,gamma)proxG(u,gamma) - u;
proxGi = @(x,gamma,i)x + (1/b)*LiS{i}( proxG_Id( Li{i}(x), b*gamma ) );

F = @(x) (1/2)*sum( sum( (Phi(x) - y).^2 ) );
E = @(x,lambda) F(x) + lambda * ( Gi{1}(x) + Gi{2}(x) + Gi{3}(x) + Gi{4}(x) );

Phis = @(f)K(M(f));

nablaF = @(x)Phis( Phi(x) - y );

beta = 1;

n = 4;

gamma = 1.8*beta;

lambda = 1e-4;

exo2()

%% Insert your code here.

x = zeros(N,N);

z = zeros(N,N,n);

exo3()

%% Insert your code here.

lambdaList = logspace( -4, -2, 10 );

exo4()

%% Insert your code here.

clf
imageplot(recov)
title( sprintf( '\\lambda=%.1e; SNR=%.2fdB', bestLambda, SNRmax ) );

n = 5;

Jmin = 4;
Psi = @(x)perform_wavortho_transf(x,Jmin,+1);

Psis = @(x)perform_wavortho_transf(x,Jmin,-1);

l1norm = @(wx)abs(wx);
proxl1 = @(wx,gamma) max(0,1 - gamma./l1norm(wx) ).*wx;
proxG5 = @(x,gamma) Psis( proxl1(Psi(x),gamma) );

exo5()

%% Insert your code here.
