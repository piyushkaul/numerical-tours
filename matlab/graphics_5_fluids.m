
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/graphics_5_fluids')

n = 128;
options.bound = 'per';
V = perform_blurring(randn(n,n,2), 40, options);

myplot = @(V)plot_vf( V(1:6:n,1:6:n, :) );

clf;
myplot(V);

normalize = @(V)V ./ repmat( max(1e-9,sqrt(sum3(V.^2, 3))) , [1 1 2]);

clf;
myplot(normalize(V));

[Y,X] = meshgrid(0:n-1,0:n-1);
mu = sin(X*pi()/n).^2; mu = -4*( mu+mu' );
mu(1) = 1;

A = @(V)real( ifft2( fft2( div(V, options) ) ./ mu ) );

ProjI = @(V)V - grad(A(V), options);

U = ProjI(V);
clf;
myplot(U);

clf;
myplot(V-U);

name = 'lena';
f = crop( load_image(name, 2*n), n);

U = normalize(ProjI(V));

periodic = @(P)cat(3, mod(P(:,:,1)-1,n)+1, mod(P(:,:,2)-1,n)+1 );

extend1 = @(f)[f f(:,1)];
extend = @(f)extend1(extend1(f)')';

myinterp = @(P1,f1,Pi)interp2( P1(:,:,2), P1(:,:,1),f1, Pi(:,:,2), Pi(:,:,1) );

[Y,X] = meshgrid(1:n,1:n); 
P = cat(3, X,Y);
[Y1,X1] = meshgrid(1:n+1,1:n+1);
P1 = cat(3, X1,Y1);

W = @(f,U)myinterp( P1, extend(f), periodic( P - U ) );

rho = 2;
clf;
imageplot(W(f,rho*U));

exo1()

%% Insert your code here.

exo2()

%% Insert your code here.

nu = 1/10;

mu = 2*nu;

Wt = @(V,U)cat(3, W(V(:,:,1),U), W(V(:,:,2),U) );

tau = .5;

V = normalize(ProjI(V));

g = f;

g = W (g,tau*U);
V = Wt(V,tau*U);

s1 = [2:n 1]; s2 = [n 1:n-1];
Delta = @(g)1/4 *( g(s1,:,:) + g(s2,:,:) + g(:,s1,:) + g(:,s2,:) ) - g;

V = V + tau*nu*Delta(V);
g = g + tau*mu*Delta(g);

V = ProjI(V);

exo3()

%% Insert your code here.
