
warning off
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('toolbox_graph')
addpath('toolbox_wavelet_meshes')
addpath('solutions/fastmarching_4_mesh')

name = 'elephant-50kv';
[vertex,faces] = read_mesh(name);
nvert = size(vertex,2);

nstart = 15;
pstarts = floor(rand(nstart,1)*nvert)+1;
options.start_points = pstarts;

clear options;
options.end_points = [];

options.W = ones(nvert,1);

options.nb_iter_max = Inf;
[D,S,Q] = perform_fast_marching_mesh(vertex, faces, pstarts, options);

clf;
plot_fast_marching_mesh(vertex,faces, D, [], options);

[Qexact,DQ, voronoi_edges] = compute_voronoi_mesh(vertex, faces, pstarts, options);
options.voronoi_edges = voronoi_edges;
plot_fast_marching_mesh(vertex,faces, D, [], options);

exo1()

%% Insert your code here.

nend = 40;
pend = floor(rand(nend,1)*nvert)+1;

vring = compute_vertex_ring(faces);

exo2()

%% Insert your code here.

options.method = 'continuous';
paths = compute_geodesic_mesh(D, vertex, faces, pend, options);

plot_fast_marching_mesh(vertex,faces, Q, paths, options);

clear options;
name = 'fandisk';
[vertex,faces] = read_mesh(name);
options.name = name;
nvert = size(vertex,2);

clf;
plot_mesh(vertex,faces, options);

options.verb = 0;
[Umin,Umax,Cmin,Cmax] = compute_curvature(vertex,faces,options);

C = abs(Cmin)+abs(Cmax);
C = min(C,.1);

options.face_vertex_color = rescale(C);
clf;
plot_mesh(vertex,faces,options);
colormap jet(256);

epsilon = .5;
W = rescale(-min(C,0.1), .1,1);

options.face_vertex_color = rescale(W);
clf;
plot_mesh(vertex,faces,options);
colormap jet(256);

pstarts = [2564; 16103; 15840];
options.start_points = pstarts;

options.W = W;
options.nb_iter_max = Inf;
[D,S,Q] = perform_fast_marching_mesh(vertex, faces, pstarts, options);

options.colorfx = 'equalize';
clf;
plot_fast_marching_mesh(vertex,faces, D, [], options);

exo3()

%% Insert your code here.

exo4()

%% Insert your code here.

exo5()

%% Insert your code here.

clear options;
options.base_mesh = 'ico';
options.relaxation = 1;
options.keep_subdivision = 0;
[vertex,faces] = compute_semiregular_sphere(7,options);
nvert = size(vertex,2);

name = 'earth';
f = load_spherical_function(name, vertex, options);
options.name = name;

pstarts = [2844; 5777];
options.start_points = pstarts;

clf;
plot_fast_marching_mesh(vertex,faces, f, [], options);
colormap gray(256);

g = load_spherical_function('earth-grad', vertex, options);

clf;
plot_fast_marching_mesh(vertex,faces, g, [], options);
colormap gray(256);

W = rescale(-min(g,10),0.01,1);

clf;
plot_fast_marching_mesh(vertex,faces, W, [], options);
colormap gray(256);

exo6()

%% Insert your code here.

exo7()

%% Insert your code here.
