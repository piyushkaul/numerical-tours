
warning on
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/denoisingwav_2_wavelet_2d')
warning off

n = 256;
name = 'hibiscus';
f0 = rescale( load_image(name,n) );
if using_matlab()
    f0 = rescale( sum(f0,3) );
end

clf;
imageplot(f0);

sigma = .08;

f = f0 + sigma*randn(size(f0));

clf;
imageplot(clamp(f), strcat(['Noisy, SNR=' num2str(snr(f0,f),3)]));

T = 1;

alpha = linspace(-3,3,1000);
clf;
plot(alpha, alpha.*(abs(alpha)>T));
axis tight;

options.ti = 0;
Jmin = 4;

a = perform_wavelet_transf(f,Jmin,+1,options);

clf;
plot_wavelet(a,Jmin);

T = 3*sigma;

aT = a .* (abs(a)>T);

clf;
plot_wavelet(aT,Jmin);

fHard = perform_wavelet_transf(aT,Jmin,-1,options);

clf;
imageplot(clamp(fHard), strcat(['Hard denoising, SNR=' num2str(snr(f0,fHard),3)]));

T = 1;
alpha = linspace(-3,3,1000);
alphaT = max(1-T./abs(alpha), 0).*alpha;
clf;
plot(alpha, alphaT);
axis tight;

T = 3/2*sigma;

aT = perform_thresholding(a,T,'soft');

aT(1:2^Jmin,1:2^Jmin) = a(1:2^Jmin,1:2^Jmin);

fSoft = perform_wavelet_transf(aT,Jmin,-1,options);

clf;
imageplot(clamp(fSoft), strcat(['Soft denoising, SNR=' num2str(snr(f0,fSoft),3)]));

exo1()

%% Insert your code here.

m = 4;

[dY,dX] = meshgrid(0:m-1,0:m-1);
delta = [dX(:) dY(:)]';

fTI = zeros(n,n);

i = 1;

fS = circshift(f,delta(:,i));

a = perform_wavelet_transf(fS,Jmin,1,options);
aT = perform_thresholding(a,T,'hard');
fS = perform_wavelet_transf(aT,Jmin,-1,options);

fS = circshift(fS,-delta(:,i));

fTI = (i-1)/i*fTI + 1/i*fS;

exo2()

%% Insert your code here.

exo3()

%% Insert your code here.

extend_stack_size(4);

options.ti = 1;
a = perform_wavelet_transf(f0,Jmin,+1,options);

clf;
i = 0;
for j=1:2
    for k=1:3
        i = i+1;
        imageplot(a(:,:,i+1), strcat(['Scale=' num2str(j) ' Orientation=' num2str(k)]), 2,3,i );
    end
end

options.ti = 1;
a = perform_wavelet_transf(f,Jmin,+1,options);

T = 3.5*sigma;
aT = perform_thresholding(a,T,'hard');

J = size(a,3)-5;
clf;
imageplot(a(:,:,J), 'Noisy coefficients');

clf;
imageplot(aT(:,:,J), 'Thresholded coefficients');

fTI = perform_wavelet_transf(aT,Jmin,-1,options);

clf;
imageplot(clamp(fTI), strcat(['Hard invariant, SNR=' num2str(snr(f0,fTI),3)]));

exo4()

%% Insert your code here.
