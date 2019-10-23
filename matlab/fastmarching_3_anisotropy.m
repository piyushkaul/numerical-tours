
warning off
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('toolbox_graph')
addpath('solutions/fastmarching_3_anisotropy')
warning on

name = 'fingerprint';
n = 150;
M = rescale(load_image(name,n));

options.order = 2;
G = grad(M,options);

T = zeros(n,n,2,2);
T(:,:,1,1) = G(:,:,2).^2;
T(:,:,2,2) = G(:,:,1).^2;
T(:,:,1,2) = -G(:,:,1).*G(:,:,2);
T(:,:,2,1) = -G(:,:,1).*G(:,:,2);

sigma = 12;
T = perform_blurring(T,sigma);

[e1,e2,l1,l2] = perform_tensor_decomp(T);

clf;
plot_vf(e1(1:10:n,1:10:n,:),M);
colormap(gray(256));

anisotropy = .1;

H = perform_tensor_recomp(e1,e2, ones(n),ones(n)*1/anisotropy );

pstart = [n n]/4;

hx = 1/n; hy = 1/n;
[D, dUx, dUy, Vor, L] = fm2dAniso([hx;hy], H, pstart);

clf;
subplot(1,2,1);
imageplot(M, 'Image');
subplot(1,2,2);
hold on;
imageplot(convert_distance_color(D), 'Geodesic distance');
hh = plot(pstart(2),pstart(1), 'r.');
set(hh, 'MarkerSize',15);
axis('ij');
colormap(gray(256));

exo1()

%% Insert your code here.

anisotropy = .02;
H = perform_tensor_recomp(e1,e2, ones(n),ones(n)*1/anisotropy );

vertex = [1;1];

[D, dUx, dUy, Vor, L] = fm2dAniso([hx;hy], H, vertex);

[tmp,i] = max(D(:));
[x,y] = ind2sub([n n],i); 
vertex(:,end+1) = [x;y];

clf;
subplot(1,2,1);
hold on;
imageplot(M, 'Image'); axis ij;
plot(vertex(2,1), vertex(1,1), 'r.');
plot(vertex(2,2), vertex(1,2), 'b.');
subplot(1,2,2);
hold on;
imageplot( convert_distance_color(D), 'Distance'); axis ij;
plot(vertex(2,1), vertex(1,1), 'r.');
plot(vertex(2,2), vertex(1,2), 'b.');
colormap gray(256);

[D1, dUx, dUy, Vor, L] = fm2dAniso([hx;hy], H, vertex);

clf;
imageplot( D, 'Old distance', 1,2,1 );
imageplot( D1, 'New distance', 1,2,2 );
colormap jet(256);

D = D1;

exo2()

%% Insert your code here.

n = 256;
name = 'peppers-bw';
M = rescale(load_image(name, n));

exo3()

%% Insert your code here.

exo4()

%% Insert your code here.

exo5()

%% Insert your code here.

exo6()

%% Insert your code here.
