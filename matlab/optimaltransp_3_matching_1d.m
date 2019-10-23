
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/optimaltransp_3_matching_1d')

n = 256;
f = rescale( load_image('lena', n) );

clf;
imageplot(f);

Q = 50;

[h,t] = hist(f(:), Q);

clf;
bar(t,h*Q/n^2); axis('tight');

exo1()

%% Insert your code here.

g = rescale( mean(load_image('fingerprint', n),3) );

clf;
imageplot(g);

exo2()

%% Insert your code here.

[~,sigmaf] = sort(f(:));
[~,sigmag] = sort(g(:));

sigmafi = [];
sigmafi(sigmaf) = 1:n^2;

sigma = sigmag(sigmafi);

f1 = reshape(g(sigma), [n n]);
imageplot(f1)

clf;
imageplot(f, 'f', 1,2,1);
imageplot(f1, '\pi_g(f)',  1,2,2);

ft = @(t)reshape( t*f1 + (1-t)*f, [n n]);

clf;
imageplot(ft(1/2));

exo3()

%% Insert your code here.


