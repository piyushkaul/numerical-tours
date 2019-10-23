
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/graphics_3_synthesis_diffusion')

n = 256;

M1 = perform_blurring(randn(n,n), 15);

x = linspace(0,1, n*n);
M2 = perform_hist_eq(M1,x);

clf;
imageplot(M1, 'Original', 1,2,1);
imageplot(M2, 'Equalized', 1,2,2);

exo1()

%% Insert your code here.

M = rescale( load_image('lena', n) );
Gr = grad(M);
tv = sum(sum( sqrt(sum(Gr.^2,3)), 2 ), 1);
disp(strcat(['Total variation of M = ' num2str(tv) '.']));

exo2()

%% Insert your code here.

exo3()

%% Insert your code here.
