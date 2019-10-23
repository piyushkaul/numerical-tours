
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/multidim_3_multispectral')

extend_stack_size(10);

name = 'unclebens';
options.nbdims = 3;
M = read_bin(name, options);
n = 256;
M = rescale( crop(M,[n n size(M,3)]) );

n = size(M,1);

p = size(M,3);

clf;
imageplot(M(:,:,10), 'R', 1,3,1);
imageplot(M(:,:,15), 'G', 1,3,2);
imageplot(M(:,:,20), 'B', 1,3,3);

rgbsel = [10 15 20];
clf;
imageplot(M(:,:,rgbsel), 'RGB');

pos = [30 50];

v = M(pos(1),pos(2),:); v = v(:);
clf;
plot( v, '.-');
axis('tight');

U = reshape( M, [n*n p] )';
U = dct(U);
U = reshape(U', [n n p]);

clf;
subplot(2,1,1);
v = M(pos(1),pos(2),:); v = v(:);
plot(v); axis('tight');
title('Spetral content');
subplot(2,1,2);
v = U(pos(1),pos(2),:); v = v(:);
plot(v); axis('tight');
title('DCT tranform');

ilist = [1 4 8 16];
clf;
for i=1:length(ilist);
    imageplot(U(:,:,ilist(i)), ['DCT freq ' num2str(ilist(i))], 2,2,i);
end

Jmin = 3;
UW = U;
for i=1:p
    UW(:,:,i) = perform_wavelet_transf(U(:,:,i), Jmin, +1);
end

clf;
subplot(1,2,1);
plot_wavelet(UW(:,:,1), Jmin);
subplot(1,2,2);
plot_wavelet(UW(:,:,10), Jmin);

m = round(.01*n*n*p);

UWT = perform_thresholding(UW, m, 'largest');

exo1()

%% Insert your code here.

exo2()

%% Insert your code here.

exo3()

%% Insert your code here.
