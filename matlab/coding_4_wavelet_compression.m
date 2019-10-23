
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/coding_4_wavelet_compression')

v = linspace(-1,1, 2048);

T = .1;

vI = floor(abs(v/T)).*sign(v);

vQ = sign(vI) .* (abs(vI)+.5) * T;

clf;
subplot(1,2,1);
plot(v, vI);
axis('tight');
title(strcat(['Quantized integer values, T=' num2str(T)]));
subplot(1,2,2);
hold('on');
plot(v, vQ); 
plot(v, v, 'r--');
axis('equal'); axis('tight');
title('De-quantized real values');

n = 256;
M = rescale( load_image('lena', n) );

Jmin = 4;
MW = perform_wavelet_transf(M,Jmin, +1);

exo1()

%% Insert your code here.

exo2()

%% Insert your code here.

exo3()

%% Insert your code here.

MWI = floor(abs(MW/T)).*sign(MW);
MWQ = sign(MWI) .* (abs(MWI)+.5) * T;

a = max(abs(MWI(:))); 
t = -a:a;
h = hist(MWI(:), t); h = h/sum(h);

t0 = 0:1/T0;
MI = floor(abs(M/T0)); % quantized pixel values
h0 = hist(MI(:), t0); h0 = h0/sum(h0);

clf;
subplot(2,1,1);
bar(t0,h0); axis('tight');
title('Pixels');
subplot(2,1,2);
bar(t,h); axis([-5 5 0 max(h)])
title('Wavelets (zoom)');

exo4()

%% Insert your code here.

exo5()

%% Insert your code here.

T = .1;
MWI = floor(abs(MW/T)).*sign(MW);

MWH = MWI(1:n/2,n/2+1:n);
MWV = MWI(n/2+1:n,1:n/2);
MWD = MWI(n/2+1:n,n/2+1:n);

clf;
imageplot(MWH,'Horizontal',1,3,1);
imageplot(MWV,'Vertical',1,3,2);
imageplot(MWD,'Diagonal',1,3,3);

exo6()

%% Insert your code here.

exo7()

%% Insert your code here.

exo8()

%% Insert your code here.
