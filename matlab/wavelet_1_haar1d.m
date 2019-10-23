
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/wavelet_1_haar1d')

N = 512;

name = 'piece-regular';
f = rescale( load_signal(name, N) );

J = log2(N)-1;

haar = @(a)[   a(1:2:length(a)) + a(2:2:length(a));
                    a(1:2:length(a)) - a(2:2:length(a)) ]/sqrt(2);

clf;
subplot(2,1,1);
plot(f); axis('tight'); title('Signal');
subplot(2,1,2);
plot(haar(f)); axis('tight'); title('Transformed');

j = J;

fw = f;

fw(1:2^(j+1)) = haar(fw(1:2^(j+1)));

s = 400; t = 40; 
clf;
subplot(2,1,1);
plot(f,'.-'); axis([s-t s+t 0 1]); title('Signal (zoom)');
subplot(2,1,2);
plot(fw(1:2^j),'.-'); axis([(s-t)/2 (s+t)/2 min(fw(1:2^j)) max(fw(1:2^j))]); title('Averages (zoom)');

exo1()

%% Insert your code here.

fprintf('Should be 0: %.3f.\n', (norm(f)-norm(fw))/norm(f));

clf; 
plot_wavelet(fw);
axis([1 N -2 2]);

ihaar = @(a,d)assign( zeros(2*length(a),1), ...
        [1:2:2*length(a), 2:2:2*length(a)], [a+d; a-d]/sqrt(2) );

f1 = fw;
j = 0;

f1(1:2^(j+1)) = ihaar(f1(1:2^j), f1(2^j+1:2^(j+1)));

exo2()

%% Insert your code here.

fprintf('Should be 0: %.3f.\n', (norm(f-f1))/norm(f));

j = J-1;

fw1 = fw; fw1(2^j+1:end) = 0;

exo3()

%% Insert your code here.

S = @(x,T) x .* (abs(x)>T);

T = .5;

fwT = S(fw,T);

clf;
subplot(2,1,1);
plot_wavelet(fw); axis('tight'); title('Original coefficients');
subplot(2,1,2);
plot_wavelet(fwT); axis('tight'); title('Thresholded coefficients');

exo4()

%% Insert your code here.

exo5()

%% Insert your code here.
