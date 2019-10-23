
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/optim_1_gradient_descent')

eta = 10;

f = @(x)( x(1)^2 + eta*x(2)^2 ) /2;

t = linspace(-.7,.7,101);
[u,v] = meshgrid(t,t);
F = ( u.^2 + eta*v.^2 )/2 ;

clf; hold on;
imagesc(t,t,F); colormap jet(256);
contour(t,t,F, 20, 'k');
axis off; axis equal;

Gradf = @(x)[x(1); eta*x(2)];

tau = 1.8/eta;

exo1()

%% Insert your code here.

clf; hold on;
imagesc(t,t,F); colormap jet(256);
contour(t,t,F, 20, 'k');
h = plot(X(1,:), X(2,:), 'k.-');
set(h, 'LineWidth', 2);
set(h, 'MarkerSize', 15);
axis off; axis equal;

exo2()

%% Insert your code here.

n = 256;
x0 = rescale( load_image('lena',n) );

clf;
imageplot(x0);

grad = @(x)cat(3, x-x([end 1:end-1],:), x-x(:,[end 1:end-1]));

v = grad(x0);

clf;
imageplot(v(:,:,1), 'd/dx', 1,2,1);
imageplot(v(:,:,2), 'd/dy', 1,2,2);

clf; 
imageplot(v);

clf; 
imageplot( sqrt( sum3(v.^2,3) ) );

div = @(v)v([2:end 1],:,1)-v(:,:,1) + v(:,[2:end 1],2)-v(:,:,2);

delta = @(x)div(grad(x));

clf; 
imageplot(delta(x0));

dotp = @(a,b)sum(a(:).*b(:));
fprintf('Should be 0: %.3i\n', dotp(grad(x0), grad(x0)) + dotp(delta(x0),x0) );

sigma = .1;
y = x0 + randn(n)*sigma;

clf;
imageplot(clamp(y));

lambda = .3/5;

epsilon = 1e-3;

NormEps = @(u,epsilon)sqrt(epsilon^2 + sum(u.^2,3));
J = @(x,epsilon)sum(sum(NormEps(grad(x),epsilon)));

f = @(y,x,epsilon)1/2*norm(x-y,'fro')^2 + lambda*J(x,epsilon);

Normalize = @(u,epsilon)u./repmat(NormEps(u,epsilon), [1 1 2]);
GradJ = @(x,epsilon)-div( Normalize(grad(x),epsilon) );

Gradf = @(y,x,epsilon)x-y+lambda*GradJ(x,epsilon);

tau = 1.8/( 1 + lambda*8/epsilon );
tau = tau*4;

exo3()

%% Insert your code here.

clf;
imageplot(clamp(x));

n = 64;
name = 'square';
x0 = load_image(name,n);

a = 4;
Lambda = ones(n);
Lambda(end/2-a:end/2+a,:) = 0;

Phi  = @(x)x.*Lambda;
PhiS = @(x)Phi(x);

y = Phi(x0);

clf;
imageplot(x0, 'Original', 1,2,1);
imageplot(y, 'Damaged', 1,2,2);

ProjH = @(x,y) x + PhiS( y - Phi(x) );

exo4()

%% Insert your code here.

exo5()

%% Insert your code here.
