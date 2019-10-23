
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/pde_4_wavelet_compression')

N = 1024;

f = load_signal('piece-polynomial', N);

clf; 
plot(f); axis tight;

phi = @(x)exp( -abs(x).^2 / 2 );

c = 1/2;

alpha = .05; 
beta = 1e-3;

sigma = @(x)alpha*abs(x-c)+beta;

normalize = @(K)K ./ repmat( sum(K,2), [1 N] );

t = (0:N-1)/N;
[Y,X] = meshgrid(t,t);
K = normalize( phi( (X-Y) ./ sigma(X)  ) );

I = round( linspace(1,N,17) );
clf;
plot(K(:,I)); axis tight;
axis([1 N -.005 .1]);

clf;
imageplot(K.^.1);

g = K*f;

clf;
subplot(2,1,1);
plot(f); axis('tight');
title('f');
subplot(2,1,2);
plot(g); axis('tight');
title('g');

t0 = 0;
ntrials = 50;
for i=1:ntrials
    tic; g = K*f; t0 = t0 + toc();
end
t0 = t0 / ntrials;

options.filter = 'haar';
Jmin = 1;

options.separable = 0;
Psi  = @(f)perform_wavelet_transf(f, Jmin, +1, options);
PsiS = @(f)perform_wavelet_transf(f, Jmin, -1, options);

options.separable = 1;
Theta = @(K)perform_wavelet_transf(K, Jmin, +1, options);

clf;
plot_wavelet(Psi(f), Jmin);

S = Theta(K);

clf;
plot_wavelet(S, Jmin, options);

Gamma = @(S,tau)S .* ( abs(S)>tau );

subs = @(u,m)u(m);
tau = @(m)subs( sort(abs(S(:)), 'descend'), m);

m = 4*N;

S1 = sparse(Gamma(S,tau(m)));

clf; 
plot_wavelet(zeros(N), Jmin, options); hold on;
spy(S1);

exo1()

%% Insert your code here.

exo2()

%% Insert your code here.
