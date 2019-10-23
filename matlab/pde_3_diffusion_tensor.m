
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/pde_3_diffusion_tensor')

n = 256;

name = 'hibiscus';
f = load_image(name,n);
f = rescale( sum(f,3) );

clf;
imageplot(f);

cconv = @(f,h)real(ifft2(fft2(f).*repmat(fft2(h),[1 1 size(f,3)])));

t = [0:n/2 -n/2+1:-1];
[X2,X1] = meshgrid(t,t);
normalize = @(h)h/sum(h(:));
h = @(sigma)normalize( exp( -(X1.^2+X2.^2)/(2*sigma^2) ) );

blur = @(f,sigma)cconv(f,h(sigma));

options.order = 2;
nabla = @(f)grad(f,options);

tensorize = @(u)cat(3, u(:,:,1).^2, u(:,:,2).^2, u(:,:,1).*u(:,:,2));

rotate = @(T)cat(3, T(:,:,2), T(:,:,1), -T(:,:,3));

T = @(f,sigma)blur( tensorize( nabla(f) ), sigma);

options.sub = 8;
clf; sigma = .1;
plot_tensor_field(rotate(T(f,sigma)), f, options);
title(['\sigma=' num2str(sigma)]);

clf; sigma = 4;
plot_tensor_field(rotate(T(f,sigma)), f, options);
title(['\sigma=' num2str(sigma)]);

clf; sigma = 10;
plot_tensor_field(rotate(T(f,sigma)), f, options);
title(['\sigma=' num2str(sigma)]);

delta = @(S)(S(:,:,1)-S(:,:,2)).^2 + 4*S(:,:,3).^2;
eigenval = @(S)deal( ...
    (S(:,:,1)+S(:,:,2)+sqrt(delta(S)))/2,  ...
    (S(:,:,1)+S(:,:,2)-sqrt(delta(S)))/2 );

normalize = @(u)u./repmat(sqrt(sum(u.^2,3)), [1 1 2]);
eig1 = @(S)normalize( cat(3,2*S(:,:,3), S(:,:,2)-S(:,:,1)+sqrt(delta(S)) ) );

ortho = @(u)deal(u, cat(3,-u(:,:,2), u(:,:,1)));
eigbasis = @(S)ortho(eig1(S));

sigma = 2;
S = T(f,sigma);
[lambda1,lambda2] = eigenval(S);
[e1,e2] = eigbasis(S);

recompose = @(lambda1,lambda2,e1,e2)repmat(lambda1,[1 1 3]).*tensorize(e1) + repmat(lambda2,[1 1 3]).*tensorize(e2);

mynorm = @(x)norm(x(:));
S1 = recompose(lambda1,lambda2,e1,e2);
fprintf('Should be 0: %.3f\n', mynorm(S-S1));

E = sqrt(lambda1+lambda2);
A = (lambda1-lambda2)./(lambda1+lambda2);

clf;
imageplot({E A}, {'E', 'A'});

m = 4;
Cm = 3.31488;

phi = @(s,lambda)1-exp( -Cm./(s/lambda).^m );

s = linspace(0,3,1024)';
clf;
plot(s, [phi(s,1) s.*phi(s,1)], 'LineWidth', 2); 
legend('\phi(s)', 's \phi(s)');

lambda = 1e-3;

sigma = 2;

S = T(f,sigma);
[lambda1,lambda2] = eigenval(S);
[e1,e2] = eigbasis(S);

S = recompose(phi(lambda1,lambda),ones(n),e1,e2);

Mult = @(S,u)cat(3, S(:,:,1).*u(:,:,1) + S(:,:,3).*u(:,:,2), ...
                          S(:,:,3).*u(:,:,1) + S(:,:,2).*u(:,:,2) );

tau = .05;

f1 = f;

f1 = f1 + tau * div( Mult(S, nabla(f1) ) );

exo1()

%% Insert your code here.
