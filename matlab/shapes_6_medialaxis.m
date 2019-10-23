
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('toolbox_graph')
addpath('solutions/shapes_6_medialaxis')

n = 200;
W = load_image('mountain', n);
W = rescale(W,.25,1);

pstart = [[20;20] [120;100] [180;30] [60;160]];
nbound = size(pstart,2);

ms = 20;
clf; hold on;
imageplot(W);
h = plot(pstart(2,:), pstart(1,:), '.r'); set(h, 'MarkerSize', ms);

[D,S,Q] = perform_fast_marching(W, pstart);

clf; hold on;
imageplot(convert_distance_color(D, W));
h = plot(pstart(2,:), pstart(1,:), '.r'); set(h, 'MarkerSize', ms);

clf; hold on;
imageplot(Q);
h = plot(pstart(2,:), pstart(1,:), '.r'); set(h, 'MarkerSize', ms);
colormap jet(256);

G = grad(Q);

G(G<-nbound/2) = G(G<-nbound/2) + nbound;
G(G>nbound/2) = G(G>nbound/2) - nbound;

G = sqrt(sum(G.^2,3));

B = 1 - (G>.1);

clf; hold on;
imageplot(B);
h = plot(pstart(2,:), pstart(1,:), '.r'); set(h, 'MarkerSize', ms);

n = 200;
name = 'chicken';
M = load_image(name,n);
M = perform_blurring(M,5);
M = double( rescale( M )>.5 );
if M(1)==1 
    M = 1-M;
end

pstart = compute_shape_boundary(M);
nbound = size(pstart,2);

lw = 2;
clf; hold on;
imageplot(-M);
h = plot(pstart(2,:), pstart(1,:), 'r'); set(h, 'LineWidth', lw); axis ij;

W = ones(n);
L = zeros(n)-Inf; L(M==1) = +Inf;

options.constraint_map = L;
[D,S,Q] = perform_fast_marching(W, pstart, options);
D(M==0) = Inf;

clf;
hold on;
display_shape_function(D);
h = plot(pstart(2,:), pstart(1,:), 'r'); set(h, 'LineWidth', lw); axis ij;

clf;
hold on;
display_shape_function(Q);
h = plot(pstart(2,:), pstart(1,:), 'r'); set(h, 'LineWidth', lw); axis ij;

exo1()

%% Insert your code here.

exo2()

%% Insert your code here.
