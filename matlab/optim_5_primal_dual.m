
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/optim_5_primal_dual')

name = 'lena';
n = 256;
f0 = load_image(name);
f0 = rescale(crop(f0,n));

clf;
imageplot(f0);

rho = .8;
Lambda = rand(n,n)>rho;

Phi = @(f)f.*Lambda;

y = Phi(f0);

clf;
imageplot(y);

K  = @(f)grad(f);
KS = @(u)-div(u);

Amplitude = @(u)sqrt(sum(u.^2,3));
F = @(u)sum(sum(Amplitude(u)));

ProxF = @(u,lambda)max(0,1-lambda./repmat(Amplitude(u), [1 1 2])).*u;

t = -linspace(-2,2, 201);
[Y,X] = meshgrid(t,t);
U = cat(3,Y,X);
V = ProxF(U,1);
% 3D display
clf;
surf(V(:,:,1)); 
colormap jet(256);
view(150,40);
axis('tight');
camlight;
shading interp;

ProxFS = @(y,sigma)y-sigma*ProxF(y/sigma,1/sigma);

V = ProxFS(U,1);
% display
clf;
surf(V(:,:,1));
colormap jet(256);
view(150,40);
axis('tight');
camlight; shading interp;

ProxG = @(f,tau)f + Phi(y - Phi(f));

L = 8;
sigma = 10;
tau = .9/(L*sigma);
theta = 1;

f = y;
g = K(y)*0;
f1 = f;

fold = f;
g = ProxFS( g+sigma*K(f1), sigma);
f = ProxG(  f-tau*KS(g), tau);
f1 = f + theta * (f-fold);

exo1()

%% Insert your code here.

clf;
imageplot(f);

exo2()

%% Insert your code here.

n = 64;
name = 'square';
f0 = load_image(name,n);

a = 4;
Lambda = ones(n);
Lambda(end/2-a:end/2+a,:) = 0;
Phi = @(f)f.*Lambda;

clf;
imageplot(f0, 'Original', 1,2,1);
imageplot(Phi(f0), 'Damaged', 1,2,2);

exo3()

%% Insert your code here.


