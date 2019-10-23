
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('toolbox_graph')
addpath('solutions/meshdeform_3_flattening')

name = 'nefertiti'; 
options.name = name;
[vertex,faces] = read_mesh(name);
n = size(vertex,2);

clf;
plot_mesh(vertex,faces, options);
shading faceted;

options.symmetrize = 1;
options.normalize = 0;
L = compute_mesh_laplacian(vertex,faces,'conformal',options);

[U,S] = eig(full(L)); S = diag(S);
[S,I] = sort(S,'ascend'); U = U(:,I);

vertexF = U(:,2:3)';

icenter = 88;
irotate = 154;
vertexF = vertexF - repmat(vertexF(:,icenter), [1 n]);
theta = -pi/2+atan2(vertexF(2,irotate),vertexF(1,irotate));
vertexF = [vertexF(1,:)*cos(theta)+vertexF(2,:)*sin(theta); ...
           -vertexF(1,:)*sin(theta)+vertexF(2,:)*cos(theta)];

clf;
plot_mesh(vertexF,faces);

exo1()

%% Insert your code here.

D = zeros(n);
for i=1:n
    D(:,i) = perform_fast_marching_mesh(vertex,faces,i);
end

D = (D+D')/2;

J = eye(n) - ones(n)/n;
W = -J*(D.^2)*J;

[U,S] = eig(W);
S = diag(S);
[S,I] = sort(S,'descend'); U = U(:,I);

clf;
hh = plot(S(1:30), '.-'); axis('tight');
set(hh, 'LineWidth', 2);

vertexF = U(:,1:2)' .* repmat(sqrt(S(1:2)), [1 n]);

vertexF = vertexF - repmat(vertexF(:,icenter), [1 n]);
theta = -pi/2+atan2(vertexF(2,irotate),vertexF(1,irotate));
vertexF = [vertexF(1,:)*cos(theta)+vertexF(2,:)*sin(theta); ...
           -vertexF(1,:)*sin(theta)+vertexF(2,:)*cos(theta)];

clf;
plot_mesh(vertexF,faces,options);

exo2()

%% Insert your code here.

clf;
plot(stress(2:end), '.-');
axis('tight');

exo3()

%% Insert your code here.
