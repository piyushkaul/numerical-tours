
warning off
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/denoisingwav_6_curvelets')
warning on

n = 256;
name = 'lena';
M = rescale(load_image(name, n));

options.null = 0;
options.finest = 1;
options.nbscales = 4;
options.nbangles_coarse = 16;
options.is_real = 1;
options.n = n;

MW = perform_curvelet_transform(M, options);

clf;
plot_curvelet(MW, options);

T = .2;
MWT = perform_thresholding(MW, T, 'hard');

M1 = perform_curvelet_transform(MWT, options);

clf;
imageplot(M, 'Original', 1,2,1);
imageplot(clamp(M1), 'Approximated', 1,2,2);

name = 'lena';
n = 128;
M0 = rescale(crop(load_image(name),n, [108 200]));
options.n = n;

sigma = .05;
M = M0 + sigma*randn(n);

exo1()

%% Insert your code here.

clf;
imageplot(clamp(M), 'Noisy', 1,2,1);
imageplot(clamp(Mcurv), ['Denoised, SNR=' num2str(snr(M0,Mcurv),3) 'dB'], 1,2,2);

exo2()

%% Insert your code here.

clf;
imageplot(clamp(M), 'Noisy', 1,2,1);
imageplot(clamp(Mcurv), ['Denoised, SNR=' num2str(snr(M0,Mcurv),3)], 1,2,2);

exo3()

%% Insert your code here.

exo4()

%% Insert your code here.
