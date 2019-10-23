
warning off
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('toolbox_graph')
addpath('solutions/fastmarching_1_2d')
warning on

clear options;
n = 300;
name = 'road2';
f = rescale( load_image(name, n) );

clf;
imageplot(f);

x0 = [14;161];
x1 = [293;148];

epsilon = 1e-2;
W = epsilon + abs(f-f(x0(1),x0(2)));

clf;
imageplot(W);

options.nb_iter_max = Inf;
options.end_points = x1;

[D,S] = perform_fast_marching(1./W, x0, options);

clf;
hold on;
imageplot( convert_distance_color(D,f) );
h = plot(x0(2),x0(1), '.r'); set(h, 'MarkerSize', 25);
h = plot(x1(2),x1(1), '.b'); set(h, 'MarkerSize', 25);

exo1()

%% Insert your code here.

options.nb_iter_max = Inf;
options.end_points = [];
[D,S] = perform_fast_marching(1./W, x0, options);

clf;
imageplot(D);
colormap jet(256);

options.order = 2;
G0 = grad(D, options);

G = G0 ./ repmat( sqrt( sum(G0.^2, 3) ), [1 1 2]);

clf;
imageplot(G);
colormap jet(256);

tau = .8;

gamma = x1;

Geval = @(G,x)[interp2(1:n,1:n,G(:,:,1),x(2),x(1)); ...
             interp2(1:n,1:n,G(:,:,2),x(2),x(1)) ];

g = Geval(G, gamma(:,end));

gamma(:,end+1) = gamma(:,end) - tau*g;

exo2()

%% Insert your code here.

clf; hold on;
imageplot(f);
h = plot(gamma(2,:),gamma(1,:), '.b'); set(h, 'LineWidth', 2);
h = plot(x0(2),x0(1), '.r'); set(h, 'MarkerSize', 25);
h = plot(x1(2),x1(1), '.b'); set(h, 'MarkerSize', 25);
axis ij;

clf; hold on;
imageplot(D); colormap jet(256);
h = plot(gamma(2,:),gamma(1,:), '.b'); set(h, 'LineWidth', 2);
h = plot(x0(2),x0(1), '.r'); set(h, 'MarkerSize', 25);
h = plot(x1(2),x1(1), '.b'); set(h, 'MarkerSize', 25);
axis ij;

exo3()

%% Insert your code here.

exo4()

%% Insert your code here.

n = 256;
name = 'cortex';
f = rescale( sum(load_image(name,n),3) );

clf;
imageplot(f);

G = grad(f,options);
G = sqrt( sum(G.^2,3) );

sigma = 3;
Gh = perform_blurring(G,sigma);

clf;
imageplot(Gh);

epsilon = 0.01;
W = 1./( epsilon + Gh );

clf;
imageplot(W);

x0 = [ [136;53] [123;205]];

options.nb_iter_max = Inf;
options.end_points = [];
[D,S,Q] = perform_fast_marching(1./W, x0, options);

clf; hold on;
imageplot( perform_hist_eq(D,'linear') );
h = plot(x0(2,:),x0(1,:), '.r'); set(h, 'MarkerSize', 25);
colormap jet(256);

clf; hold on;
A = zeros(n,n,3); A(:,:,1) = rescale(Q); A(:,:,3) = f;
imageplot(A);
h = plot(x0(2,:),x0(1,:), '.g'); set(h, 'MarkerSize', 25);

exo5()

%% Insert your code here.

exo6()

%% Insert your code here.

n = 256;
name = 'vessels';
f = rescale(load_image(name, n));

clf;
imageplot(f);

sigma = 20;
f1 = perform_blurring(f,sigma) - f;

clf;
imageplot(f1);

c = max(f1(:));
epsilon = 1e-2;
W = epsilon + abs(f1-c);

clf,
imageplot(W);

x0 = [142;226];

exo7()

%% Insert your code here.

exo8()

%% Insert your code here.

x0 = [[143;249] [174;9]];

exo9()

%% Insert your code here.
