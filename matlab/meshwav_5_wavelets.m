
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('toolbox_graph')
addpath('toolbox_wavelet_meshes')
addpath('solutions/meshwav_5_wavelets')

options.base_mesh = 'ico';
options.relaxation = 1;
options.keep_subdivision = 1;
J = 6;
[vertex,face] = compute_semiregular_sphere(J,options);

options.use_color = 1;
options.rho = .3;
options.color = 'rescale';
options.use_elevation = 0;

f = load_spherical_function('earth', vertex{end}, options);

clf;
plot_spherical_function(vertex,face,f, options);
colormap gray(256);

fw = perform_wavelet_mesh_transform(vertex,face, f, +1, options);

r = .1;
fwT = perform_thresholding( fw, round(r*length(fw)), 'largest' );

f1 = perform_wavelet_mesh_transform(vertex,face, fwT, -1, options);

clf;
subplot(1,2,1);
plot_spherical_function(vertex,face,f, options);
title('Original function');
subplot(1,2,2);
plot_spherical_function(vertex,face,f1, options);
title('Approximated function');
colormap gray(256);

exo1()

%% Insert your code here.

exo2()

%% Insert your code here.

exo3()

%% Insert your code here.

name = 'bunny';
M = read_gim([name '-sph.gim']);
n = size(M,1);

clf;
imageplot(M);

clf;
plot_geometry_image(M, 1,1);
view(20,88);

options.order = 2;
u = zeros(n,n,3); v = zeros(n,n,3);
for i=1:3
    [u(:,:,i),v(:,:,i)] = grad(M(:,:,i), options);
end

v = cat(3, u(:,:,2).*v(:,:,3)-u(:,:,3).*v(:,:,2), ...
    u(:,:,3).*v(:,:,1)-u(:,:,1).*v(:,:,3), ...
    u(:,:,1).*v(:,:,2)-u(:,:,2).*v(:,:,1) );

L = [1 2 -1]; L = reshape(L/norm(L), [1 1 3]);
A1 = max( sum( v .* repmat(L, [n n]), 3 ), 0 );
L = [-1 -2 -1]; L = reshape(L/norm(L), [1 1 3]);
A2 = max( sum( v .* repmat(L, [n n]), 3 ), 0 );

clf;
imageplot(A1, '', 1,2,1);
imageplot(A2, '', 1,2,2);

J = 6;
[vertex,face,vertex0] = compute_semiregular_gim(M,J,options);

options.func = 'mesh';
options.name = name;
options.use_elevation = 0;
options.use_color = 0;

selj = J-3:J;
clf;
for j=1:length(selj)
    subplot(2,2,j);
    plot_mesh(vertex{selj(j)},face{selj(j)}, options);
    shading('faceted'); lighting('flat'); axis tight;
    % title(['Subdivision level ' num2str(selj(j))]);
end
colormap gray(256);

f = vertex{end}';

fw = perform_wavelet_mesh_transform(vertex,face, f, +1, options);

r = .1;
fwT = perform_thresholding( fw, round(r*length(fw)), 'largest' );

f1 = perform_wavelet_mesh_transform(vertex,face, fwT, -1, options);

clf;
subplot(1,2,1);
plot_mesh(f,face{end},options); shading('interp'); axis('tight');
title('Original surface');
subplot(1,2,2);
plot_mesh(f1,face{end},options); shading('interp'); axis('tight');
title('Wavelet approximation');
