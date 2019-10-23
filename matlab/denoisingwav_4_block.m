
warning off
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/denoisingwav_4_block')
warning on

n = 256;

name = 'boat';
f0 = rescale( load_image(name,n) );

clf; imageplot(f0);

sigma = .08;

f = f0 + sigma*randn(size(f0));

clf; imageplot(clamp(f));

Jmin = 4;
options.ti = 0;

wav  = @(f)perform_wavelet_transf(f,Jmin,+1,options);
iwav = @(fw)perform_wavelet_transf(fw,Jmin,-1,options);

clf;
plot_wavelet(wav(f),Jmin);

psi= @(s,T)max3(1-T^2 ./ max(abs(s).^2,1e-9),0);

t = linspace(-3,3,1024);
clf; hold on;
plot(t,t.*psi(t,1)); 
plot(t,t, 'r--'); axis equal;

theta = @(x,T)psi(x,T).*x;
ThreshWav = @(f,T)iwav(theta(wav(f),T));

T = 1.5*sigma;
clf;
imageplot(clamp( ThreshWav(f,T) ));

exo1()

%% Insert your code here.

clf;
imageplot(clamp(fThresh), strcat(['SNR=' num2str(snr(f0,fThresh),3)]));

w = 4;

[dX,dY,X,Y] = ndgrid(0:w-1,0:w-1,1:w:n-w+1,1:w:n-w+1);
I = X+dX + (Y+dY-1)*n;

block = @(x)reshape(x(I(:)),size(I));

iblock = @(H)assign(zeros(n), I, H);

mynorm = @(x)norm(x(:));
fprintf('Should be 0: %.3f\n', mynorm(f - iblock(block(f))) );

repm = @(v)repmat( max3(v,1e-15), [w w]);
energy = @(H)repm( sqrt( mean(mean(abs(H).^2,1),2) ) );

Thresh = @(H,T)psi(energy(H),T).*H;
ThreshBlock = @(x,T)iblock( Thresh(block(x),T) );

exo2()

%% Insert your code here.

T = 1.25*sigma;
clf;
plot_wavelet( ThreshBlock(wav(f),T), Jmin);

ThreshWav = @(f,T)iwav(ThreshBlock(wav(f),T));

clf;
imageplot(clamp( ThreshWav(f,T) ));

exo3()

%% Insert your code here.

clf;
imageplot(clamp(fBlock), strcat(['SNR=' num2str(snr(f0,fBlock),3)]));

options.ti = 1;
wav  = @(f)perform_wavelet_transf(f,Jmin,+1,options);
iwav = @(fw)perform_wavelet_transf(fw,Jmin,-1,options);

fw = wav(f);

[dX,dY,X,Y,J] = ndgrid(0:w-1,0:w-1,1:w:n-w+1,1:w:n-w+1, 1:size(fw,3));
I = X+dX + (Y+dY-1)*n + (J-1)*n^2;

block = @(x)reshape(x(I(:)),size(I));
iblock = @(H)assign(zeros(size(fw)), I, H);

repm = @(v)repmat( max3(v,1e-15), [w w]);
energy = @(H)repm( sqrt( mean(mean(abs(H).^2,1),2) ) );

Thresh = @(H,T)psi(energy(H),T).*H;
ThreshBlock = @(x,T)iblock( Thresh(block(x),T) );

ThreshWav = @(f,T)iwav(ThreshBlock(wav(f),T));

clf;
T = 1.25*sigma;
imageplot(clamp( ThreshWav(f,T) ));

exo4()

%% Insert your code here.

clf;
imageplot(clamp(fTI), strcat(['SNR=' num2str(snr(f0,fTI),3)]));
