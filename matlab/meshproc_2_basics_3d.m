
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('toolbox_graph')
addpath('solutions/meshproc_2_basics_3d')

name = 'elephant-50kv';
options.name = name; % useful for displaying
[vertex,face] = read_mesh(name);

clf;
plot_mesh(vertex, face, options);
shading('interp');

clf;
for i=1:4
    subplot(2,2,i);
    plot_mesh(vertex, face);
    shading('faceted');
    zoom(1.8^(i+1));
end

options.face_vertex_color = vertex(1,:)';
clf;
plot_mesh(vertex, face, options);
colormap(jet(256));
camlight;

options.face_vertex_color = cos(50*vertex(2,:)');
clf;
plot_mesh(vertex, face, options);
colormap(jet(256));

options.face_vertex_color = [];

noise = randn(size(vertex))*.01;
noise(:,vertex(1,:)>mean(vertex(1,:))) = 0;
vertex1 = vertex+noise;

clf;
subplot(1,2,1);
plot_mesh(vertex,face, options); axis('tight');
subplot(1,2,2);
plot_mesh(vertex1,face, options); axis('tight');

vertex1 = sign(vertex) .* (abs(vertex)/max(abs(vertex(:)))).^1.8;

clf;
subplot(1,2,1);
plot_mesh(vertex,face, options); axis('tight');
subplot(1,2,2);
plot_mesh(vertex1,face, options); axis('tight');

name = 'mushroom';
options.name = name;
[vertex,face] = read_mesh(name);

[normal,normalf] = compute_normal(vertex,face);

clf; 
options.normal = normal;
plot_mesh(vertex,face,options); shading('interp'); 
axis('tight');
options.normal = [];

name = 'elephant-50kv';
options.name = name;
[vertex,face] = read_mesh(name);

[normal,normalf] = compute_normal(vertex,face);
vertex1 = vertex + .02*normal;
vertex2 = vertex + .04*normal;

options.face_vertex_color =  [];
clf;
subplot(1,2,1);
plot_mesh(vertex1,face,options); shading('interp'); axis('tight');
subplot(1,2,2);
plot_mesh(vertex2,face,options); shading('interp'); axis('tight');

laplacian_type = 'distance';

options.symmetrize = 0;
options.normalize = 1;
L = compute_mesh_laplacian(vertex,face,laplacian_type,options);

options.symmetrize = 0;
options.normalize = 1;
L0 = compute_mesh_laplacian(vertex,face,laplacian_type,options);

v1 = L*vertex(1,:)';
v2 = L*vertex(2,:)';

vmax = median(abs(v1)*5); v1 = clamp(v1,-vmax,vmax);
vmax = median(abs(v2)*5); v2 = clamp(v2,-vmax,vmax);

clf;
subplot(1,2,1);
options.face_vertex_color = v1;
plot_mesh(vertex,face, options);
title('Laplacian of X');
subplot(1,2,2);
options.face_vertex_color = v2;
plot_mesh(vertex,face, options);
title('Laplacian of Y');
options.face_vertex_color = [];
colormap(jet(256));

name = 'elephant-50kv';
options.name = name; % useful for displaying
[vertex,face] = read_mesh(name);

options.curvature_smoothing = 10;
options.verb = 0;
[Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(vertex,face,options);

clf;
subplot(1,2,1);
options.face_vertex_color = perform_saturation(Cgauss,1.2);
plot_mesh(vertex,face, options); shading interp; colormap jet(256);
title('Gaussian curvature');
subplot(1,2,2);
options.face_vertex_color = perform_saturation(abs(Cmin)+abs(Cmax),1.2);
plot_mesh(vertex,face, options); shading interp; colormap jet(256);
title('Total curvature');

[vertex,face] = read_tet('toolbox_additional/hand.tet');

clf;
plot_mesh(vertex,face,options);
