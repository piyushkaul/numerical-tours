
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('toolbox_graph')
addpath('solutions/fastmarching_7_sampling_compr')

name = 'cameraman';
n = 256;
M = rescale( load_image(name, n) );

m = 400;

vertex = floor( rand(2,m-4)*(n-1) ) +1;
vertex(:,end+1:end+4) = [[1;1] [1;n] [n;n] [n;1]];

faces = compute_delaunay(vertex);

vinterp = interp2(M, vertex(2,:), vertex(1,:));

Minterp = compute_triangulation_interpolation(faces,vertex,vinterp, n);

clf;
subplot(1,2,1);
plot_triangulation(vertex,faces, M);
title('Triangulation');
subplot(1,2,2);
imageplot(clamp(Minterp), ['Interpolation, SNR=' num2str(snr(Minterp,M)) 'dB']);

vapprox = compute_orthoproj_triangulation(vertex, faces, M);
Mapprox = compute_triangulation_interpolation(faces,vertex,vapprox, n);

clf;
imageplot(clamp(Minterp), ['Interpolation, SNR=' num2str(snr(Minterp,M),3) 'dB'], 1,2,1);
imageplot(clamp(Mapprox), ['Approximation, SNR=' num2str(snr(Mapprox,M),3) 'dB'], 1,2,2);

alpha = .7;
epsilon = 1e-2;

exo1()

%% Insert your code here.

exo2()

%% Insert your code here.

vgeod = compute_orthoproj_triangulation(vertex, faces, M);
Mgeod = compute_triangulation_interpolation(faces,vertex,vgeod, n);

clf;
subplot(1,2,1);
plot_triangulation(vertex,faces, M);
subplot(1,2,2);
imageplot(clamp(Mgeod), ['SNR=' num2str(snr(Mgeod,M),3) 'dB']);

exo3()

%% Insert your code here.
