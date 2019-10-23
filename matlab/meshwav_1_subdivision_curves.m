
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('toolbox_graph')
addpath('toolbox_wavelet_meshes')
addpath('solutions/meshwav_1_subdivision_curves')

ms = 20; lw = 1.5;
myplot = @(f,c)plot(f([1:end 1]), c, 'LineWidth', lw, 'MarkerSize', ms);
myaxis = @(rho)axis([-rho 1+rho -rho 1+rho], 'off');

f0 =    [0.11 0.18 0.26 0.36 0.59 0.64 0.80 0.89 0.58 0.22 0.18 0.30 0.58 0.43 0.42]' + ...
   1i * [0.91 0.55 0.91 0.58 0.78 0.51 0.81 0.56 0.10 0.16 0.35 0.42 0.40 0.24 0.31]';

f0 = rescale(real(f0),.01,.99) + 1i * rescale(imag(f0),.01,.99);

clf; myplot(f0, 'k.-'); 
myaxis(0);

subdivide = @(f,h)cconvol( upsampling(f), h);

h = [1 4 6 4 1];
h = 2*h/sum(h(:));

f = f0;

f = subdivide(f,h);

clf; hold on;
myplot(f, 'k.-');
myplot(f0, 'r.--');
myaxis(0);

exo1()

%% Insert your code here.

exo2()

%% Insert your code here.

exo3()

%% Insert your code here.

hcc = @(w)[1 w w 1]/(1+w);

exo4()

%% Insert your code here.

exo5()

%% Insert your code here.

h4pt = @(w)[-w, 0, 1/2+w, 1, 1/2+w, 0, -w];

exo6()

%% Insert your code here.

exo7()

%% Insert your code here.

exo8()

%% Insert your code here.

H = {   [0.5000 0.5000], ...
        [-0.0625, 0.5625, 0.5625, -0.0625], ...
        [0.0117, -0.0977, 0.5859, 0.5859, -0.0977, 0.0117], ...
        [-0.0024, 0.0239, -0.1196, 0.5981, 0.5981, -0.1196, 0.0239, -0.0024] };    
hdd = @(k)assign(assign(zeros(4*k-1,1),1:2:4*k-1,H{k}), 2*k, 1);

exo9()

%% Insert your code here.

options.bound = 'per';
n = 1024*2; 
sigma = n/8;
F = perform_blurring(randn(n,1),sigma,options) + 1i*perform_blurring(randn(n,1),sigma,options);
F = rescale(real(F),.01,.99) + 1i * rescale(imag(F),.01,.99);

clf; myplot(F, 'k');
myaxis(0);

h = [-1, 0, 9, 1, 9, 0, -1]/16; 
h((end+1)/2)=1;

exo10()

%% Insert your code here.

dist = @(f,g)abs( repmat(f, [1 length(g)]) - repmat(transpose(g), [length(f) 1]) );

hausdorff = @(f,g)sqrt( mean(min(dist(f,g)).^2) );
hausdorff = @(f,g)hausdorff(f,g) + hausdorff(g,f);

exo11()

%% Insert your code here.

% exo12()

%% Insert your code here.
