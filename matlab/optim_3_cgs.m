
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/optim_3_cgs')

n = 500;

A = randn(n);
A = A*A' + .1*eye(n);

b = randn(n,1);

dotp = @(a,b)sum(a(:).*b(:));

exo1()

%% Insert your code here.

n = 256;
g = rescale( load_image('lena',n) );

clf;
imageplot(g);

s = [n 1:n-1];
grad = @(f)cat(3, f-f(s,:), f-f(:,s));

v = grad(g);

clf;
imageplot(v(:,:,1), 'd/dx', 1,2,1);
imageplot(v(:,:,2), 'd/dy', 1,2,2);

clf; 
imageplot(v);

clf; 
imageplot( sqrt( sum3(v.^2,3) ) );

t = [2:n 1];
div = @(v)v(t,:,1)-v(:,:,1) + v(:,t,2)-v(:,:,2);

delta = @(f)div(grad(f));

clf; 
imageplot(delta(g));

dotp = @(a,b)sum(a(:).*b(:));
fprintf('Should be 0: %.3i\n', dotp(grad(g), grad(g)) + dotp(delta(g),g) );

M = rand(n)>.7;
w = 30;
M(end/4-w:end/4+w,end/4-w:end/4+w) = 0;

Phi = @(x)M.*x;

y = Phi(g);

clf;
imageplot(y);

lambda = .01;

A = @(x)Phi(x) - lambda*delta(x);

b = y;

exo2()

%% Insert your code here.

clf;
imageplot(clamp(x));

A = @(z)cat(3, -delta(z(:,:,1)) + Phi(z(:,:,2)), Phi(z(:,:,1)) );

b = cat(3, zeros(n), y);

exo3()

%% Insert your code here.

clf;
imageplot(x);
