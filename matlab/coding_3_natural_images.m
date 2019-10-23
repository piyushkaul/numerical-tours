
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/coding_3_natural_images')

n = 256;
M1 = rescale( load_image('boat', n) );
M2 = rescale( load_image('lena', n) );

clf;
subplot(2,2,1);
imageplot(M1);
subplot(2,2,2);
[h,t] = hist(M1(:), 60);
bar(t, h/sum(h));
axis('square');
subplot(2,2,3);
imageplot(M2);
subplot(2,2,4);
[h,t] = hist(M2(:), 60);
bar(t, h/sum(h));
axis('square');

exo1()

%% Insert your code here.

n = 256*2;
M = rescale( load_image('lena', n) );

Jmin = 4;
MW = perform_wavelet_transf(M,Jmin, +1);

MW1 = MW(1:n/2,n/2+1:n);

v = max(abs(MW1(:)));
k = 20;
t = linspace(-v,v,2*k+1);
h = hist(MW1(:), t);

clf;
subplot(1,2,1);
imageplot(MW1);
subplot(1,2,2);
bar(t, h/sum(h)); axis('tight');
axis('square');

T = .03;
MW1q = floor(abs(MW1/T)).*sign(MW1);
nj = size(MW1,1);

t = 2; % you can try with other values
C = reshape(MW1q([ones(1,t) 1:nj-t],:),size(MW1));

options.normalize = 1;
[H,x,xc] = compute_conditional_histogram(MW1q,C, options);

q = 8; % width for display
H = H((end+1)/2-q:(end+1)/2+q,(end+1)/2-q:(end+1)/2+q);
clf;
imagesc(-q:q,-q:q,max(log(H), -5)); axis image;
colormap gray(256);

options.normalize = 0;
[H,x,xc] = compute_conditional_histogram(MW1q,C, options);
H = H((end+1)/2-q:(end+1)/2+q,(end+1)/2-q:(end+1)/2+q);

clf;
contour(-q:q,-q:q,max(log(H), -5), 20, 'k'); axis image;
colormap gray(256);

C = sign( reshape(MW1q([1 1:nj-1],:),size(MW1)) );

[H,x,xc] = compute_conditional_histogram(MW1q,C);

clf; plot(x,H, '.-');
axis([-10 10 0 max(H(:))]);
legend('sign=-1', 'sign=0', 'sign=+1');
set_graphic_sizes([], 20);

ent_total = compute_entropy(MW1q);

h0 = compute_histogram(C);
H(H==0) = 1e-9; % avoid numerical problems
ent_partial = sum( -log2(H).*H );
ent_cond = sum( ent_partial.*h0' );

disp(['Global coding: ' num2str(ent_total,3), ' bpp.']);
disp(['Conditional coding: ' num2str(ent_cond,3), ' bpp.']);
