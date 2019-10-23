
warning off
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('toolbox_graph')
addpath('toolbox_graph\mex')
addpath('solutions\fastmarching_1_2d')
warning on

clear options;
options.order = 2;

n = 321;
name = 'road2';
%f = rescale( load_image(name, n) );

vars= load('vars');
T=permute(vars.gmetric, [3,4,1,2]);
metric1d = zeros(size(T,1), size(T,2));
for ix=1:size(T,1)
    for jx=1:size(T,2)
        mtx = T(ix,jx,:,:);
        mtx = squeeze(mtx);
        dist = [1 1]*mtx*[1;1];
        metric1d(ix,jx) = dist;
    end 
end 

minmetric = min(metric1d(:));
metric1d =  metric1d;% - minmetric;
figure;surf(metric1d);

T=double(T);
%T=T(75+1:75+150,75+1:75+150,:,:);
f=rescale(vars.z, n);%(160-75+1:160+75,160-75+1:160+75)
G = grad(f,options);
orig_f=f;


clf;
imageplot(f);

x0 = [154;139];
x1 = [320;320];

epsilon = 1e-2;
%W = epsilon + abs(f-f(x0(1),x0(2)));
W=T;

figure;
clf;
imageplot(metric1d);

options.nb_iter_max = Inf;
options.end_points = x1;
options.nb_iter_max = 10;

W2 = ones(300,300,2,2);
W2 = T;
W = metric1d;
%[D,S] = perform_fast_marching(1./W, x0, options);
[D,S] = perform_fast_marching(W, x0, options);
%[D, dUx, dUy, Vor, L] = fm2dAniso([1;1], W2, [1;1]);

figure;
clf;
hold on;
imageplot( convert_distance_color(D,f) );
h = plot(x0(2),x0(1), '.r'); set(h, 'MarkerSize', 25);
h = plot(x1(2),x1(1), '.b'); set(h, 'MarkerSize', 25);
colorbar

exo1()

%% Insert your code here.

options.nb_iter_max = Inf;
options.end_points = [];
[D,S] = perform_fast_marching(W, x0, options);

figure;
clf;
imageplot(D);
colormap jet(256);

options.order = 2;
G0 = grad(D, options);

G = G0 ./ repmat( sqrt( sum(G0.^2, 3) ), [1 1 2]);

figure;
clf;
imageplot(G);
colormap jet(256);

tau = .8;

gamma = x1;

%Geval = @(G,x)[interp2(1:n,1:n,G(:,:,1,1),x(2),x(1)) interp2(1:n,1:n,G(:,:,1,2),x(2),x(1)); ...
%             interp2(1:n,1:n,G(:,:,2,1),x(2),x(1)) interp2(1:n,1:n,G(:,:,2,2),x(2),x(1))];
         
Geval = @(G,x)[interp2(1:n,1:n,G(:,:,1),x(2),x(1)); ...
             interp2(1:n,1:n,G(:,:,2),x(2),x(1)) ];         

g = Geval(G, gamma(:,end));

gamma(:,end+1) = gamma(:,end) - tau*g;

exo2()

%% Insert your code here.

clf; hold on;
imageplot(f);
h = plot(gamma(2,:),gamma(1,:), '.b'); set(h, 'LineWidth', 2);
h = plot(x0(2),x0(1), '.r'); set(h, 'MarkerSize', 25);
h = plot(x1(2),x1(1), '.b'); set(h, 'MarkerSize', 25);
axis ij;

clf; hold on;
imageplot(D); colormap jet(256);
h = plot(gamma(2,:),gamma(1,:), '.b'); set(h, 'LineWidth', 2);
h = plot(x0(2),x0(1), '.r'); set(h, 'MarkerSize', 25);
h = plot(x1(2),x1(1), '.b'); set(h, 'MarkerSize', 25);
axis ij;

%exe3_nn()

%% Insert your code here.

exo4()

%% Insert your code here.

n = 321;
name = 'cortex';
f = W;%rescale( sum(load_image(name,n),3) );

clf;
imageplot(f);

G = grad(f,options);
G = sqrt( sum(G.^2,3) );

sigma = 3;
Gh = perform_blurring(G,sigma);

clf;
imageplot(Gh);

epsilon = 0.01;
%W = 1./( epsilon + Gh );

clf;
imageplot(W);

x0 = [ [136;53] [123;205]];

options.nb_iter_max = Inf;
options.end_points = [];
[D,S,Q] = perform_fast_marching(W, x0, options);

clf; hold on;
imageplot( perform_hist_eq(D,'linear') );
h = plot(x0(2,:),x0(1,:), '.r'); set(h, 'MarkerSize', 25);
colormap jet(256);

clf; hold on;
A = zeros(n,n,3); A(:,:,1) = rescale(Q); A(:,:,3) = f;
imageplot(A);
h = plot(x0(2,:),x0(1,:), '.g'); set(h, 'MarkerSize', 25);

exo5()

%% Insert your code here.

exo6()

%% Insert your code here.

n = 256;
name = 'vessels';
f = rescale(load_image(name, n));

clf;
imageplot(f);

sigma = 20;
f1 = perform_blurring(f,sigma) - f;

clf;
imageplot(f1);

c = max(f1(:));
epsilon = 1e-2;
W = epsilon + abs(f1-c);

clf,
imageplot(W);

x0 = [142;226];

exo7()

%% Insert your code here.

exo8()

%% Insert your code here.

x0 = [[2;321] [174;165]];
W = metric1d;
f=orig_f;
n=321;
exo9()

%% Insert your code here.
