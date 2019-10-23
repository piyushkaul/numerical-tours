
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('toolbox_graph')
addpath('solutions/fastmarching_8_segmentation')

n = 256;
name = 'cortex';
M = rescale( sum(load_image(name,n),3) );

clf;
imageplot(M);

pstart = [154;175];

W = abs(M-M(pstart(1),pstart(2)));
W = rescale( max(W,0.03), 0.01,1).^2;

clear options;
options.nb_iter_max = Inf;
options.end_points = [];
[D,S,Q] = perform_fast_marching(1./W, pstart, options);

exo1()

%% Insert your code here.

mu = 2;
d = sqrt( sum( grad(perform_blurring(M,mu)).^2, 3) ); 
d = perform_blurring(d,mu);

W = rescale( min(d,0.15), 0.01,1).^2;

clf;
imageplot(W);

pstart = [[30;30] [139;86] [158;170] [128;134] [124;122]];

options.nb_iter_max = Inf;
options.end_points = [];
[D,S,Q] = perform_fast_marching(1./W, pstart, options);

clf;
imageplot(Q);
colormap(jet(256));

exo2()

%% Insert your code here.
