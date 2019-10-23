
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/graphics_2_synthesis_wavelets')

n = 512;
name = 'texture';
M = load_image(name, n);
M = rescale( sum(M,3) );

extend_stack_size(4);

options.ti = 1;
Jmin = 4;
MW = perform_wavelet_transf(M(:,:,1), Jmin, +1, options);

M1 = perform_hist_eq(randn(n,n), M);

clf;
imageplot(M, 'Exemplar', 1,2,1);
imageplot(M1, 'Initial noise', 1,2,2);

MW1 = perform_wavelet_transf(M1, Jmin, +1, options);

for i=1:size(MW,3)
    MW1(:,:,i) = perform_hist_eq(MW1(:,:,i), MW(:,:,i));    
end

M1 = perform_wavelet_transf(MW1, Jmin, -1, options);

clf;
imageplot(M, 'Exemplar', 1,2,1);
imageplot(M1, 'Initial synthesis', 1,2,2);

exo1()

%% Insert your code here.

n = 512;
M = rescale( load_image('texture', n) );

exo2()

%% Insert your code here.

exo3()

%% Insert your code here.

M1 = randn(n,n,3);

[U,R] = qr(randn(3));

d = reshape(M,[n^2 3])*U;
d1 = reshape(M1,[n^2 3])*U;

for c=1:3
    d1(:,c) = perform_hist_eq(d1(:,c),d(:,c));
end

M1 = reshape(d1*U',[n n 3]);

m = M(:,:,1); m1 = M1(:,:,1);
clf;
subplot(2,1,1);
hist(m(:),50); title('Original');
subplot(2,1,2);
hist(clamp(m1(:)),50); title('Matched');

M1 = randn(n,n,3);
for i=1:3
    M1(:,:,i) = perform_hist_eq(M1(:,:,i), M(:,:,i));
end

exo4()

%% Insert your code here.

m = M(:,:,1); m1 = M1(:,:,1);
clf;
subplot(2,1,1);
hist(m(:),50); title('Original');
subplot(2,1,2);
hist(clamp(m1(:)),50); title('Matched');

clf;
imageplot(M, 'Image', 1,2,1);
imageplot(M1, 'Equalized', 1,2,2);

exo5()

%% Insert your code here.
