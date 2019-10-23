
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('toolbox_graph')
addpath('solutions/fastmarching_5_sampling_2d')

n = 256; % size of the image
sigma = n/8; % width of the bumps
[Y,X] = meshgrid(1:n,1:n);
x = n/4; y = n/4;
M = exp( -( (X-x).^2 + (Y-y).^2 )/(2*sigma^2) );
x = 3*n/4; y = 3*n/4;
M = M + exp( -( (X-x).^2 + (Y-y).^2 )/(2*sigma^2) );

W = rescale(M,1e-2,1);

m = 20; 
pstart = floor( rand(2,m)*(n-1) ) +1;

[D,Z,Q] = perform_fast_marching(1./W, pstart);

clf;
subplot(1,2,1);
hold on;
imageplot(perform_hist_eq(D,'linear'));  title('Geodesic distance'); 
plot(pstart(2,:), pstart(1,:), 'r.');
subplot(1,2,2);
hold on;
imageplot(Q); title('Voronoi'); 
plot(pstart(2,:), pstart(1,:), 'r.');
colormap jet(256);

exo1()

%% Insert your code here.

clf;
subplot(1,2,1);
hold on;
imageplot(Q, 'Voronoi'); axis ij;
plot(pstart(2,:), pstart(1,:), 'r.');
subplot(1,2,2);
plot_triangulation(pstart,faces, M);
colormap jet(256);

faces_euc = compute_delaunay(pstart);

clf;
subplot(1,2,1);
plot_triangulation(pstart,faces, M); title('Geodesic');
subplot(1,2,2);
plot_triangulation(pstart,faces_euc, M); title('Euclidean');
colormap jet(256);

W = rescale(M,3*1e-1,1);

vertex = [1;1];

[D,Z,Q] = perform_fast_marching(1./W, vertex);

[tmp,i] = max(D(:));
[x,y] = ind2sub([n n],i); 
vertex(:,end+1) = [x;y];

clf;
subplot(1,2,1);
hold on;
imageplot(W, 'Metric'); axis ij;
plot(vertex(2,1), vertex(1,1), 'r.');
plot(vertex(2,2), vertex(1,2), 'b.');
subplot(1,2,2);
hold on;
imageplot( perform_hist_eq(D, 'linear'), 'Distance'); axis ij;
plot(vertex(2,1), vertex(1,1), 'r.');
plot(vertex(2,2), vertex(1,2), 'b.');
colormap jet(256);

options.constraint_map = D;
[D1,Z,Q] = perform_fast_marching(1./W, vertex(:,end), options);

clf;
imageplot( convert_distance_color(D,M,1), 'Update distance', 1,3,1 );
imageplot( convert_distance_color(D1,M,1), 'Update distance', 1,3,2 );
imageplot( convert_distance_color(min(D,D1),M,1), 'New distance', 1,3,3 );

D = min(D,D1);

exo2()

%% Insert your code here.

exo3()

%% Insert your code here.

n = 512;
W = ones(n);

p = 40;
vertex = floor(rand(2,p)*n-1)+1;

[D,Z,Q] = perform_fast_marching(1./W, vertex);

clf; hold on;
imageplot(Q');
h = plot(vertex(1,:), vertex(2,:), 'k.');
set(h, 'MarkerSize', 15);
colormap(jet(256));

for i=1:p
    [x,y] = ind2sub(size(W), find(Q==i));
    vertex(:,i) = [mean(x);mean(y)];
end

[D,Z,Q] = perform_fast_marching(1./W, vertex);
clf; hold on;
imageplot(Q');
h = plot(vertex(1,:), vertex(2,:), 'k.');
set(h, 'MarkerSize', 15);
colormap(jet(256));

exo4()

%% Insert your code here.

n = 256;
x = linspace(-1,1,n);
[Y,X] = meshgrid(x,x);
sigma = .3;
W = exp( -(X.^2+Y.^2)/(2*sigma^2) );
W = rescale(W,.02,1);

clf;
imageplot(W);

exo5()

%% Insert your code here.
