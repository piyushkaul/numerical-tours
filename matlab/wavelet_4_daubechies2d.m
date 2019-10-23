
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/wavelet_4_daubechies2d')

p = 4;

[h,g] = compute_wavelet_filter('Daubechies',p);

disp(['h filter = [' num2str(h) ']']);
disp(['g filter = [' num2str(g) ']']);

N = 256;
f = rand(N,1);

a = subsampling( cconvol(f,h) );
d = subsampling( cconvol(f,g) );

f1 =  cconvol(upsampling(a),reverse(h)) + cconvol(upsampling(d),reverse(g));

disp(strcat((['Error |f-f1|/|f| = ' num2str(norm(f-f1)/norm(f))])));

n = 256;
name = 'hibiscus';
f = load_image(name,n);
f = rescale( sum(f,3) );

j = log2(n)-1;

fW = f;

A = fW(1:2^(j+1),1:2^(j+1));

Coarse = subsampling(cconvol(A,h,1),1);
Detail = subsampling(cconvol(A,g,1),1);

A = cat3(1, Coarse, Detail );

clf;
imageplot(f,'Original imge',1,2,1);
imageplot(A,'Vertical transform',1,2,2);

Coarse = subsampling(cconvol(A,h,2),2);
Detail = subsampling(cconvol(A,g,2),2);

A = cat3(2, Coarse, Detail );

fW(1:2^(j+1),1:2^(j+1)) = A;

clf;
imageplot(f,'Original image',1,2,1);
subplot(1,2,2);
plot_wavelet(fW,log2(n)-1); title('Transformed')

exo1()

disp(strcat(['Energy of the signal       = ' num2str(norm(f(:)).^2)]));
disp(strcat(['Energy of the coefficients = ' num2str(norm(fW(:)).^2)]));

clf;
subplot(1,2,1);
imageplot(f); title('Original');
subplot(1,2,2);
plot_wavelet(fW, 1); title('Transformed');

f1 = fW;
j = 0;

A = f1(1:2^(j+1),1:2^(j+1));

Coarse = A(1:2^j,:);
Detail = A(2^j+1:2^(j+1),:);

Coarse = cconvol(upsampling(Coarse,1),reverse(h),1);
Detail = cconvol(upsampling(Detail,1),reverse(g),1);

A = Coarse + Detail;

Coarse = A(:,1:2^j);
Detail = A(:,2^j+1:2^(j+1));

Coarse = cconvol(upsampling(Coarse,2),reverse(h),2);
Detail = cconvol(upsampling(Detail,2),reverse(g),2);

A = Coarse + Detail;

f1(1:2^(j+1),1:2^(j+1)) = A;

exo2()

disp(strcat((['Error |f-f1|/|f| = ' num2str(norm(f(:)-f1(:))/norm(f(:)))])));

eta = 4;
fWLin = zeros(n,n);
fWLin(1:n/eta,1:n/eta) = fW(1:n/eta,1:n/eta);

exo3()

T = .2;

fWT = fW .* (abs(fW)>T);

clf;
subplot(1,2,1);
plot_wavelet(fW); axis('tight'); title('Original coefficients');
subplot(1,2,2);
plot_wavelet(fWT); axis('tight'); title('Thresholded coefficients');

exo4()

exo5()

%% Insert your code here.

exo6()

clf;
subplot(1,2,1);
opt.separable = 0;
plot_wavelet(fW,1,opt);
title('Isotropic wavelets');
subplot(1,2,2);
opt.separable = 1;
plot_wavelet(fWSep,1,opt);
title('Separable wavelets');

exo7()

disp(strcat((['Error |f-f1|/|f| = ' num2str(norm(f(:)-f1(:))/norm(f(:)))])));

exo8()

exo9()

exo10()
