
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/inverse_6_image_separation')

n = 256;
name = 'boat';
f = rescale( crop(load_image(name),n) );

clf;
imageplot(f);

u = linspace(-5,5)';
clf;
subplot(2,1,1); hold('on');
plot(u, abs(u), 'b');
plot(u, sqrt(.5^2+u.^2), 'r');
title('\epsilon=1/2'); axis('square');
subplot(2,1,2); hold('on');
plot(u, abs(u), 'b');
plot(u, sqrt(1^2+u.^2), 'r');
title('\epsilon=1'); axis('square');

epsilon = 1e-2;

J = @(u)sum(sum( sqrt( epsilon^2 + sum3(grad(u).^2,3) ) ));
disp(['J(f) = ' num2str(J(f),3)]);

lambda = .2;

tau = 1.9 / ( 1 + lambda * 8 / epsilon);

u = f;

GradJ0 = @(Gr)-div( Gr ./ repmat( sqrt( epsilon^2 + sum3(Gr.^2,3) ) , [1 1 2]) );
GradJ = @(u)GradJ0(grad(u));

u = u - tau*( u - f + lambda* GradJ(u) );

exo1()

%% Insert your code here.

clf;
imageplot(u);

rho = .8; % constrast factor 
eta = .2; % saturation limit
displaytexture0 = @(x)sign(x).*abs(x).^rho;
displaytexture  = @(v)displaytexture0( clamp(v,-eta,eta)/eta );

clf;
imageplot( displaytexture(f-u) );

eta = .05;
x = [0:n/2-1, -n/2:-1]/n;
[Y,X] = meshgrid(x,x);
W = 1 ./ (eta + sqrt(X.^2 + Y.^2));

imageplot(fftshift(1./W));

T = @(v)1/2*norm( W.*fft2(v)/n, 'fro' ).^2;
disp(['T(f) = ' num2str(T(f), 3) ] );

GradT = @(f)real(ifft2(W.^2.*fft2(f)));

imageplot(GradT(f));

GradTInv = @(f)real(ifft2(fft2(f)./W.^2));

imageplot(GradTInv(f));

lambda = 5;

tau = 1.9 /( max(W(:))^2 + 8*lambda/epsilon );

u = f;

u = u - tau * ( GradT(u-f) + lambda*GradJ(u) );

exo2()

%% Insert your code here.

clf;
imageplot(u);

clf;
imageplot( displaytexture(f-u) );

p = [125 200];

mu = 10;

[Y,X] = meshgrid( 1:n, 1:n );
U = exp( ( -(X-p(1)).^2 - (Y-p(2)).^2 )/(2*mu^2)  );

clf;
imageplot(U.*f);

F = fft2(U.*f);
F = fftshift(F);
F(end/2-20:end/2+20,end/2-20:end/2+20) = 0;

[tmp,i] = max(abs(F(:))); 
[xm,ym] = ind2sub([n n], i);

clf; hold on;
imageplot(abs(F));
h = plot( [ym n-ym], [xm  n-xm], 'r.' );
set(h, 'MarkerSize', 20);

r = sqrt( (xm-n/2)^2 + (ym-n/2)^2 );

sigma = 10;
x = [0:n/2-1, -n/2:-1];
[Y,X] = meshgrid(x,x);
R = sqrt(X.^2+Y.^2);
W = 1 - exp( -(r-R).^2 / (2*sigma^2) );

imageplot(fftshift(W));

exo3()

%% Insert your code here.

exo4()

%% Insert your code here.

clf;
imageplot(u);

clf;
imageplot( displaytexture(f-u) );
