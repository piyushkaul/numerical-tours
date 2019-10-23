
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/optim_2_newton')

f = @(x1,x2)(1-x1).^2 + 100*(x2-x1.^2).^2;

x1 = linspace(-2,2,150);
x2 = linspace(-.5,3,150);
[X2,X1] = meshgrid(x2,x1);
F = f(X1,X2);

clf; 
surf(x2,x1, F, perform_hist_eq(F, 'linear') ); 
% shading interp;
% camlight;
% axis tight;
% colormap jet;

clf;
imageplot( perform_hist_eq(F, 'linear') );
colormap jet(256);

gradf = @(x1,x2)[2*(x1-1) + 400*x1.*(x1.^2-x2); 200*(x2-x1.^2)];
Gradf = @(x)gradf(x(1),x(2));

hessf = @(x1,x2)[2 + 400*(x1.^2-x2) + 800*x1.^2, -400*x1; ...
                -400*x1,  200];
Hessf = @(x)hessf(x(1),x(2));

exo1()

%% Insert your code here.

exo2()

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

sigma = .1;
y = g + randn(n)*sigma;

clf;
imageplot(clamp(y));

lambda = .3/5;

epsilon = 1e-3;

Amplitude = @(u)sqrt(epsilon^2 + sum(u.^2,3));
J = @(x)sum(sum(Amplitude(grad(x))));

f = @(x)1/2*norm(x-y,'fro')^2 + lambda*J(x);

Normalize = @(u)u./repmat(Amplitude(u), [1 1 2]);
GradJ = @(x)-div( Normalize(grad(x)) );

Gradf = @(x)x-y+lambda*GradJ(x);

A = @(x)sqrt(sum(grad(x).^2, 3));
delta = @(a,u)u./sqrt( epsilon^2 + repmat(a.^2,[1 1 2]) );

H = @(z,a)z - lambda*div( delta(a, grad(z)) );

flat = @(x)x(:); flatI = @(x)reshape(x,[n,n]);
tol = eps; maxit = 40;
Hinv = @(u,a,u_prev)cgs(@(z)flat(H(flatI(z),a)), flat(u),tol,maxit,[],[],flat(u_prev));

x = y;

d = zeros(n); % replace thie line by the previous iterate d value
[d,~] = Hinv(Gradf(x), A(x), d);
d = flatI(d);

x = x - d;

exo3()

%% Insert your code here.

clf;
imageplot(clamp(x));

exo4()

%% Insert your code here.

exo5()

%% Insert your code here.
