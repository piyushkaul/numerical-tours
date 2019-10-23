
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('toolbox_graph')
addpath('solutions/shapes_3_bendinginv_landmarks')

name = 'elephant-50kv';
options.name = name;
[vertex,faces] = read_mesh(name);
nverts = size(vertex,2);

clf;
plot_mesh(vertex,faces, options);

landmarks = 23057;
Dland = [];

[Dland(:,end+1),S,Q] = perform_fast_marching_mesh(vertex, faces, landmarks(end));

[tmp,landmarks(end+1)] = max( min(Dland,[],2) );

[Dland(:,end+1),S,Q] = perform_fast_marching_mesh(vertex, faces, landmarks(end));

clf;
options.start_points = landmarks;
plot_fast_marching_mesh(vertex,faces, min(Dland,[],2) , [], options);

exo1()

%% Insert your code here.

D = Dland(landmarks,:);
D = (D+D')/2;

J = eye(n) - ones(n)/n;
K = -1/2 * J*(D.^2)*J;

opt.disp = 0; 
[Xstrain, val] = eigs(K, 3, 'LR', opt);
Xstrain = Xstrain .* repmat(sqrt(diag(val))', [n 1]);
Xstrain = Xstrain';

vertex1 = zeros(nverts,3);
deltan = mean(Dland.^2,1);
for i=1:nverts
    deltax = Dland(i,:).^2;
    vertex1(i,:) = 1/2 * ( Xstrain * ( deltan-deltax )' )';
end
vertex1 = vertex1';

clf;
plot_mesh(vertex1,faces,options);

exo2()

%% Insert your code here.
