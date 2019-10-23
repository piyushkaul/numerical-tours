
warning off
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('toolbox_graph')
addpath('solutions/fastmarching_2_3d')
warning on

name = 'vessels';
options.nbdims = 3;
M = read_bin(name, options);
M = rescale(M);

n = size(M,1);

clf;
imageplot(M(:,:,50), 'X/Y slice', 1, 3, 1);
imageplot(squeeze(M(:,50,:)), 'X/Z slice', 1, 3, 2);
imageplot(squeeze(M(50,:,:)), 'Y/Z slice', 1, 3, 3);

slices = round(linspace(10,n-10,4));
clf;
for i=1:length(slices)
    s = slices(i);
    imageplot( M(:,:,s), strcat(['Z=' num2str(s)]), 2,2,i );
end

clf;
h = vol3d('cdata',M,'texture','2D');
view(3); axis off;

colormap bone(256);
% set up an alpha map
options.sigma = .08; % control the width of the non-transparent region
options.center = .4; % here a value in [0,1]
a = compute_alpha_map('gaussian', options); % you can plot(a) to see the alphamap
% refresh the rendering
% vol3d(h);

exo1()

%% Insert your code here.

sel = 1:2:n;
clf;
isosurface( M(sel,sel,sel), .5);
axis('off');

delta = 5;
start_point = [107;15;delta];

W = abs( M - M(start_point(1),start_point(2),start_point(3)) );

W = rescale(W,1e-2,1);

options.nb_iter_max = Inf;
[D,S] = perform_fast_marching(1./W, start_point, options);

clf;
imageplot(D(:,:,delta), '', 1,2,1);
imageplot(D(:,:,n-delta), '', 1,2,2);
colormap(jet(256));

exo2()

%% Insert your code here.

options.method = 'discrete';
minpath = compute_geodesic(D,end_point,options);

Dend = D(end_point(1),end_point(2),end_point(3));
D1 = double( D<=Dend );
% clf;
% plot_fast_marching_3d(M,D1,minpath,start_point,end_point);

exo3()

%% Insert your code here.


