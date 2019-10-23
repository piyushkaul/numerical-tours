
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/denoisingadv_6_nl_means')

n = 128;

c = [100 200];
f0 = load_image('lena');
f0 = rescale( crop(f0,n, c) );

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

patch = @(f)f(X + (Y-1)*n);

P = patch(f);

clf;
for i=1:16
    x = floor( rand*(n-1)+1 );
    y = floor( rand*(n-1)+1 );
    imageplot( squeeze(P(x,y,:,:)), '', 4,4,i );
end

d = 25;

resh = @(P)reshape(P, [n*n w1*w1])';

remove_mean = @(Q)Q - repmat(mean(Q), [w1*w1 1]);

P1 = remove_mean(resh(P));
C = P1*P1';

[V,D] = eig(C); D = diag(D); 
[D,I] = sort(D, 'descend'); V = V(:,I);

clf;
plot(D); axis('tight');

clf;
for i=1:16
    imageplot( reshape(V(:,i),[w1 w1]), '', 4,4,i );
end

iresh = @(Q)reshape(Q', [n n d]);
descriptor = @(f)iresh( V(:,1:d)' * remove_mean(resh(P)) );

H = descriptor(f);

distance = @(i)sum( (H - repmat(H(i(1),i(2),:), [n n 1])).^2, 3 )/(w1*w1);

normalize = @(K)K/sum(K(:));
kernel = @(i,tau)normalize( exp( -distance(i)/(2*tau^2) ) );

tau = .05;
i = [83 72];
D = distance(i);
K = kernel(i,tau);

clf;
imageplot(D, 'D', 1,2,1);
imageplot(K, 'K', 1,2,2);

q = 14;

selection = @(i){clamp(i(1)-q:i(1)+q, 1,n), clamp(i(2)-q:i(2)+q,1,n)};

distance = @(i,sel)sum( (H(sel{1},sel{2},:) - repmat(H(i(1),i(2),:), ...
        [length(sel{1}) length(sel{2}) 1])).^2, 3 )/(w1*w1);
distance = @(i)distance(i,selection(i));
kernel = @(i,tau)normalize( exp( -distance(i)/(2*tau^2) ) );

D = distance(i);
K = kernel(i,tau);

clf;
imageplot(D, 'D', 1,2,1);
imageplot(K, 'K', 1,2,2);

NLval = @(K,sel)sum(sum(K.*f(sel{1},sel{2})));
NLval = @(i,tau)NLval( kernel(i,tau), selection(i) );

[Y,X] = meshgrid(1:n,1:n);
NLmeans = @(tau)arrayfun(@(i1,i2)NLval([i1 i2],tau), X,Y);

tau = .03;
clf; 
imageplot(NLmeans(tau));

exo1()

%% Insert your code here.

clf;
imageplot(clamp(fNL));

exo2()

%% Insert your code here.
