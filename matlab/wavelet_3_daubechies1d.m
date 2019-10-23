
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/wavelet_3_daubechies1d')

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

name = 'piece-regular';
N = 512;
f = rescale( load_signal(name, N) );

clf;
plot(f);
axis('tight');

fw = f;

j = log2(N)-1;

a1 = fw(1:2^(j+1));

a = subsampling(cconvol(a1,h));
d = subsampling(cconvol(a1,g));

fw(1:2^(j+1)) = cat(1, a, d );

clf;
subplot(2,1,1);
plot(f); axis('tight'); title('Signal');
subplot(2,1,2);
plot(fw); axis('tight'); title('Transformed');

s = 400; t = 40; 
clf;
subplot(2,1,1);
plot(f,'.-'); axis([s-t s+t 0 1]); title('Signal (zoom)');
subplot(2,1,2);
plot(a,'.-'); axis([(s-t)/2 (s+t)/2 min(a) max(a)]); title('Averages (zoom)');

exo1()

%% Insert your code here.

disp(strcat(['Energy of the signal       = ' num2str(norm(f).^2,3)]));
disp(strcat(['Energy of the coefficients = ' num2str(norm(fw).^2,3)]));

clf; 
plot_wavelet(fw);
axis([1 N -1 1]);

f1 = fw;
j = 0;

a = f1(1:2^j);
d = f1(2^j+1:2^(j+1));

a = cconvol(upsampling(a,1),reverse(h),1) + cconvol(upsampling(d,1),reverse(g),1);

f1(1:2^(j+1)) = a;

exo2()

%% Insert your code here.

disp(strcat((['Error |f-f1|/|f| = ' num2str(norm(f-f1)/norm(f))])));

T = .5;

fwT = fw .* (abs(fw)>T);

clf;
subplot(2,1,1);
plot_wavelet(fw); axis([1 N -1 1]); title('Original coefficients');
subplot(2,1,2);
plot_wavelet(fwT); axis([1 N -1 1]); title('Thresholded coefficients');

exo3()

%% Insert your code here.

exo4()

%% Insert your code here.

exo5()

%% Insert your code here.

exo6()

%% Insert your code here.
