
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/segmentation_1_edge_detection')

n = 256*2;

f0 = load_image('hibiscus',n);
f0 = rescale(sum(f0,3));

clf;
imageplot(f0);

cconv = @(f,h)real(ifft2(fft2(f).*fft2(h)));

t = [0:n/2 -n/2+1:-1];
[X2,X1] = meshgrid(t,t);
normalize = @(h)h/sum(h(:));
h = @(sigma)normalize( exp( -(X1.^2+X2.^2)/(2*sigma^2) ) );

blur = @(f,sigma)cconv(f,h(sigma));

exo1()

%% Insert your code here.

s = [n 1:n-1];
nabla = @(f)cat(3, f-f(s,:), f-f(:,s));

v = nabla(f0);

clf;
imageplot(v(:,:,1), 'd/dx', 1,2,1);
imageplot(v(:,:,2), 'd/dy', 1,2,2);

sigma = 1;
d = sqrt( sum(nabla(  blur(f0,sigma)  ).^2,3) );

clf;
imageplot(d);

exo2()

%% Insert your code here.

exo3()

%% Insert your code here.

t = [2:n 1];
div = @(v)v(t,:,1)-v(:,:,1) + v(:,t,2)-v(:,:,2);

delta = @(f)div(nabla(f));

clf; 
imageplot(delta(f0));

dotp = @(a,b)sum(a(:).*b(:));
fprintf('Should be 0: %.3i\n', dotp(nabla(f0), nabla(f0)) + dotp(delta(f0),f0) );

sigma = 4;
clf;
plot_levelset( delta(blur(f0,sigma)) ,0,f0);

exo4()

%% Insert your code here.

dx = @(f)(f(s,:)-f(t,:)) / 2;
dy = @(f)dx(f')';

s = [2:n 1]; t = [n 1:n-1];
d2x = @(f)f(s,:) + f(t,:) - 2*f;
d2y = @(f)d2x(f')';
dxy = @(f)dy(dx(f));

hessian = @(f)cat(3, d2x(f), dxy(f), d2y(f));

sigma = 6;
g = grad( blur(f0,sigma) );

h = hessian( blur(f0,sigma) );

a = h(:,:,1:2).*repmat(g(:,:,1), [1 1 2]) + ...
    h(:,:,2:3).*repmat(g(:,:,2), [1 1 2]);

clf;
plot_levelset( sum(a.*g, 3) ,0,f0);

exo5()

%% Insert your code here.
