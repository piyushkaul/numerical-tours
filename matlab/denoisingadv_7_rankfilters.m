
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/denoisingadv_7_rankfilters')

n = 256;

name = 'hibiscus';
f0 = load_image(name, n);
f0 = rescale(crop( sum(f0,3) ,n));

clf;
imageplot(f0);

sigma = .04;

f = f0 + randn(n,n)*sigma;

clf;
imageplot(clamp(f));

w = 3;
w1 = 2*w+1;

[Y,X] = meshgrid(1:n,1:n);
[dY,dX] = meshgrid(-w:w,-w:w);
dX = reshape(dX, [1 1 w1 w1]);
dY = reshape(dY, [1 1 w1 w1]);
X = repmat(X, [1 1 w1 w1]) + repmat(dX, [n n 1 1]);
Y = repmat(Y, [1 1 w1 w1]) + repmat(dY, [n n 1 1]);

X(X<1) = 2-X(X<1); Y(Y<1) = 2-Y(Y<1);
X(X>n) = 2*n-X(X>n); Y(Y>n) = 2*n-Y(Y>n);

Pi = @(f)reshape( f(X + (Y-1)*n), [n n w1*w1] );

P = Pi(f);

clf;
for i=1:16
    x = floor( rand*(n-1)+1 );
    y = floor( rand*(n-1)+1 );
    imageplot( reshape(P(x,y,:,:), w1,w1), '', 4,4,i );
end

Pmean = @(f)mean(Pi(f),3);

clf;
imageplot(Pmean(f));

p = 100;
psi = @(f)f.^(1/p);
ipsi = @(f)f.^p;
imageplot(Pmean(abs(f)) - ipsi(Pmean(psi(abs(f)))));

r = @(beta)min(ceil(beta*w1*w1)+1,w1*w1);

subsample = @(x,s)x(:,:,s);
phi = @(f,beta)subsample(sort(Pi(f), 3), r(beta));

exo1()

%% Insert your code here.

closing = @(f)phi(f,0);
clf;
imageplot(closing(f));

opening = @(f)phi(f,1);
clf;
imageplot(opening(f));

exo2()

%% Insert your code here.

exo3()

%% Insert your code here.

exo4()

%% Insert your code here.

medfilt = @(f)phi(f,1/2);

clf;
imageplot(medfilt(f));

exo5()

%% Insert your code here.

clf;
imageplot(f1);
