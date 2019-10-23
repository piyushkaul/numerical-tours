
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/multidim_2_volumetric')

extend_stack_size(10);

name = 'vessels';
options.nbdims = 3;
M = read_bin(name, options);
M = rescale(M);

n = size(M,1);

slices = round(linspace(10,n-10,4));
clf;
for i=1:length(slices)
    s = slices(i);
    imageplot( M(:,:,s), strcat(['Z=' num2str(s)]), 2,2,i );
end

sel = 1:2:n;
clf;
isosurface( M(sel,sel,sel), .5);
axis('off');

MW = M;

MW = cat3(1, (MW(1:2:n,:,:)+MW(2:2:n,:,:))/sqrt(2), (MW(1:2:n,:,:)-MW(2:2:n,:,:))/sqrt(2) );

MW = cat3(2, (MW(:,1:2:n,:)+MW(:,2:2:n,:))/sqrt(2), (MW(:,1:2:n,:)-MW(:,2:2:n,:))/sqrt(2) );

MW = cat3(3, (MW(:,:,1:2:n)+MW(:,:,2:2:n))/sqrt(2), (MW(:,:,1:2:n)-MW(:,:,2:2:n))/sqrt(2) );

clf;
imageplot(MW(:,:,30), 'Horizontal slice', 1,2,1);
imageplot(squeeze(MW(:,30,:)), 'Vertical slice', 1,2,2);

exo1()

%% Insert your code here.

m = round(.01*n^3);
MWT = perform_thresholding(MW, m, 'largest');

exo2()

%% Insert your code here.

s = 30;
clf;
imageplot( M(:,:,s), 'Original', 1,2,1 );
imageplot( clamp(M1(:,:,s)), 'Approximation', 1,2,2 );

sel = 1:2:n;
clf;
isosurface( M1(sel,sel,sel), .5);
axis('off');

sigma = .06;
Mnoisy = M + sigma*randn(n,n,n);

clf;
imageplot(Mnoisy(:,:,n/2),'X slice', 1,2,1);
imageplot(squeeze(Mnoisy(:,n/2,:)),'Y slice', 1,2,2);

x = -n/2:n/2-1;
[X,Y,Z] = ndgrid(x,x,x);

s = 2; % width
h = exp( -(X.^2 + Y.^2 + Z.^2)/(2*s^2) );
h = h/sum(h(:));

Mh = real( ifftn( fftn(Mnoisy) .* fftn(fftshift(h)) ) );

i = 40;
clf;
imageplot( Mnoisy(:,:,i), 'Noisy', 1,2,1 );
imageplot( Mh(:,:,i), 'Denoised', 1,2,2 );

sel = 1:2:n;
clf;
isosurface( Mh(sel,sel,sel), .5);
axis('off');

exo3()

%% Insert your code here.

sel = 1:2:n;
clf;
isosurface( Mblur(sel,sel,sel), .5);
axis('off');
title(['Filtering, SNR=' num2str(snr(M,Mblur),3) 'dB']);

exo4()

%% Insert your code here.

sel = 1:2:n;
clf;
isosurface( Mwav(sel,sel,sel), .5);
title(['Soft thresholding, SNR=' num2str(snr(M,Mwav),3) 'dB']);
axis('off');

w = 4;

[dX,dY,dZ] = ndgrid(0:w-1,0:w-1,0:w-1);

Mspin = zeros(n,n,n);

for i=1:w^3
    % shift the image
    MnoisyC = circshift(Mnoisy, [dX(i) dY(i) dZ(i)]);
    % denoise the image to get a result M1
    M1 = MnoisyC; % replace this line by some denoising
    % shift inverse
    M1 = circshift(M1, -[dX(i) dY(i) dZ(i)]);
    % average the result
    Mspin = Mspin*(i-1)/i + M1/i;
end

exo5()

%% Insert your code here.

sel = 1:2:n;
clf;
isosurface( Mspin(sel,sel,sel), .5);
title(['Cycle spinning, SNR=' num2str(snr(M,Mspin),3) 'dB']);
axis('off');

[h,g] = compute_wavelet_filter('Daubechies',2*4);

MW = M;
j = log2(n)-1;

A = MW(1:2^(j+1),1:2^(j+1),1:2^(j+1));
for d=1:3
    a = cat(d, subsampling(cconv(A,h,d),d), subsampling(cconv(A,g,d),d) );
end
MW(1:2^(j+1),1:2^(j+1),1:2^(j+1)) = A;

exo6()

%% Insert your code here.

MWT = perform_thresholding(MW,m,'largest');

clf;
subplot(1,2,1);
plot_wavelet(MW(:,:,n/8));
subplot(1,2,2);
plot_wavelet(MWT(:,:,n/8));

M1 = MWT;
j = 0;

A = M1(1:2^(j+1),1:2^(j+1),1:2^(j+1));
for d=1:3
    W = subselectdim(A,2^j+1:2^(j+1),d);
    A = subselectdim(A,1:2^j,d);
    A = cconv(upsampling(A,d),reverse(h),d) + cconv(upsampling(W,d),reverse(g),d);
end
M1(1:2^(j+1),1:2^(j+1),1:2^(j+1)) = A;

exo7()

%% Insert your code here.

sel = 1:2:n;
clf;
isosurface( M1(sel,sel,sel), .5);
title(['Soft thresholding, SNR=' num2str(snr(M,M1),3) 'dB']);
axis('off');

exo8()

%% Insert your code here.

sel = 1:2:n;
clf;
isosurface( Mwav(sel,sel,sel), .5);
title(['Soft thresholding, SNR=' num2str(snr(M,Mwav),3) 'dB']);
axis('off');

exo9()

%% Insert your code here.

sel = 1:2:n;
clf;
isosurface( Mspin(sel,sel,sel), .5);
title(['Cycle spinning, SNR=' num2str(snr(M,Mspin),3) 'dB']);
axis('off');
