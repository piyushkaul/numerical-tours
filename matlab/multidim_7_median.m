
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/multidim_7_median')

d = 2; % dimension
n = 1000; % number of points
X = randn(d,n);

p = 100; % number of outliers
sel = randperm(n); sel = sel(1:p); % index of outliers
X(:,sel) = rand(d,p)*50;

m = mean(X,2);

niter = 30;

med = m;
energy = [];
for i=1:niter
    % comute the distance from med to the points
    dist = sqrt( sum( (X-repmat(med,[1 n])).^2 ) );
    % compute the weight, take care of not dividing by 0
    weight = 1./max( dist, 1e-10 ); weight = weight/sum(weight);
    % compute the weighted mean
    med = sum( repmat(weight,[d 1]).*X, 2 );
    energy(end+1) = sum( dist );
end

clf;
plot(energy, '.-'); axis('tight')
set_label('Iteration', 'L1 energy');

clf;
hold('on');
plot(X(1,:), X(2,:), '.');
plot(m(1,:), m(2,:), 'k*');
plot(med(1,:), med(2,:), 'ro');
axis('tight');

name = 'flowers';
options.nbdims = 3;
n = 256;
M0 = load_image(name, n, options);
M0 = rescale(M0);

rho = .4;

mask = repmat(rand(n,n)<rho, [1 1 3]);

sigma1 = .03; sigma2 = 1;

noise = sigma1*randn(n,n,3).*(1-mask) + sigma2*rand(n,n,3).*mask;

M = M0+noise;
pnoisy = snr(M0,M);

clf;
imageplot(M0, 'Clean image', 1,2,1);
imageplot(clamp(M), strcat(['Noisy, SNR=' num2str(pnoisy)]), 1,2,2 );

k = 4;
w = 2*k+1;

exo1()

%% Insert your code here.

x = 100; y = 73;

selx = x-k:x+k; selx = mod(selx-1,n)+1;
sely = y-k:y+k; sely = mod(sely-1,n)+1;

patch = M(selx,sely,:);

X = reshape( patch, [w*w 3])';

exo2()

%% Insert your code here.

m = mean(X, 2);
clf;
hold('on');
plot3(X(1,:), X(2,:), X(3,:), '.');
plot3(m(1), m(2), m(3), '*k');
plot3(med(1), med(2), med(3), 'or');
view(3); axis('tight');

exo3()

%% Insert your code here.
