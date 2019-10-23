
warning off
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/denoisingwav_5_data_dependent')
warning on

lambda = [4 10 20];
[k,Lambda] = meshgrid(1:50, lambda);
P = Lambda.^k .* exp(-Lambda)./factorial(k);
h = plot(P'); axis('tight');
if using_matlab()
    set(h, 'LineWidth', 2);
end
legend('\lambda=2', '\lambda=10', '\lambda=20');
set_label('k', 'P(k)');

n = 256;
name = 'lena';
f0u = rescale( load_image(name,n) );

lmin = 1;
lmax = 40;
f0 = floor( rescale(f0u,lmin,lmax) );

f = poissrnd(f0);

clf;
imageplot(f0, 'Intensity map f0', 1,2,1);
imageplot(f,  'Observed noisy image f', 1,2,2);

clf;
imageplot(f0, 'Intensity map f0', 1,2,1);
imageplot(f-f0,  'f-f0', 1,2,2);

exo1()

%% Insert your code here.

Jmin = 4;
options.ti = 1;

exo2()

%% Insert your code here.

clf;
imageplot(f0, 'Original image', 1,2,1);
imageplot(clamp(fPoisson,min(f0(:)),max(f0(:))), strcat(['Denoised, SNR=' num2str(snr(f0,fPoisson),4)]), 1,2,2);

exo3()

%% Insert your code here.

exo4()

%% Insert your code here.

clf;
imageplot(clamp(fPoisson,min(f0(:)),max(f0(:))), ...
        strcat(['Un-stabilized, SNR=' num2str(snr(f0,fPoisson),4)]), 1,2,1);
imageplot(clamp(fVST,min(f0(:)),max(f0(:))), ...
        strcat(['Stabilized, SNR=' num2str(snr(f0,fVST),4)]), 1,2,2);

exo5()

%% Insert your code here.

n = 256;
name = 'boat';
f0 = rescale( load_image(name,n), 1e-2,1 );

K = 4;
sigma = 1/sqrt(K);

W = gamrnd(1/sigma^2, sigma^2,n,n);

f = f0.*W;

clf;
imageplot(f0, 'Intensity map f0', 1,2,1);
imageplot(f,  'Observed noisy image f', 1,2,2);

clf;
imageplot(f0, 'Intensity map f0', 1,2,1);
imageplot(f-f0,  'f-f0', 1,2,2);

exo6()

%% Insert your code here.

exo7()

%% Insert your code here.

clf;
imageplot(f0, 'Original image', 1,2,1);
imageplot(clamp(fMult,min(f0(:)),max(f0(:))), strcat(['Denoised, SNR=' num2str(snr(f0,fMult),4)]), 1,2,2);

a = psi(K) - log(K);

clf;
imageplot(f, 'f', 1,2,1);
imageplot(clamp(log(f)-a,-2,2), 'log(f)', 1,2,2);

clf;
subplot(2,1,1);
hist(W(:),100);
axis('tight');
subplot(2,1,2);
hist(log(W(:))-a,100);
axis('tight');

exo8()

%% Insert your code here.

clf;
imageplot(clamp(fMult,min(f0(:)),max(f0(:))), strcat(['Un-stabilized, SNR=' num2str(snr(f0,fMult),4)]), 1,2,1);
imageplot(clamp(fVST,min(f0(:)),max(f0(:))), strcat(['Stabilized, SNR=' num2str(snr(f0,fVST),4)]), 1,2,2);
