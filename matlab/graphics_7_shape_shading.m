
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('toolbox_graph')
addpath('solutions/graphics_7_shape_shading')

name = 'mozart';
f = load_image(name);
n = size(f,1);

h = n*.4;

f = rescale(f,0,h);

clf;
imageplot(f);

clf;
surf(f);
colormap(gray(256));
shading interp;
axis equal; 
view(110,45);
axis('off');
camlight;

options.order = 2; 
N0 = cat(3, -grad(f,options), ones(n));

s = sqrt( sum(N0.^2,3) );
N = N0 ./ repmat( s, [1 1 3] );

clf;
imageplot(N);

d = [0 0 1];

L = max(0, sum( N .* repmat(reshape(d,[1 1 3]), [n n 1]),3 ) );

vmin = .3;
clf;
imageplot(max(L,vmin));

exo1()

%% Insert your code here.

d = [0 0 1];
L = sum( N .* repmat(reshape(d,[1 1 3]), [n n 1]),3 );

epsilon = 1e-9;

clf;
imageplot(L>1-epsilon);

W = sqrt(1./L.^2-1);

W = max(W,epsilon);

clf;
imageplot(min(W,3));

p = [140;120];

[f1,S] = perform_fast_marching(1./W, p);

f1 = -f1*n;

clf;
hold on;
surf(f1);
h = plot3(p(2), p(1), f1(p(1),p(2)), 'r.');
set(h, 'MarkerSize', 30);
colormap(gray(256));
shading interp;
axis('equal'); 
view(110,45);
axis('off');
camlight;

exo2()

%% Insert your code here.

exo3()

%% Insert your code here.
