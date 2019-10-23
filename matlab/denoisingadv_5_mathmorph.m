
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/denoisingadv_5_mathmorph')

n = 256;
M = rescale( load_image('cortex',n) );

clf;
imageplot(M);

M = double(M>.45);

clf;
imageplot(M);

wmax = 7;
[Y,X] = meshgrid(-wmax:wmax, -wmax:wmax);
normalize = @(x)x/sum(x(:));
strel = @(w)normalize( double( X.^2+Y.^2<=w^2 ) );

exo1()

%% Insert your code here.

dillation=@(x,w)double(perform_convolution(x,strel(w))>0);
Md = dillation(M,2);

clf;
imageplot(Md);

exo2()

%% Insert your code here.

errosion=@(x,w)double( perform_convolution(x,strel(w))>=.999 );
Me = errosion(M,2);

clf;
imageplot(Me);

exo3()

%% Insert your code here.

opening = @(x,w)dillation(errosion(x,w),w);

w = 1;
Mo = opening(M,w);

clf;
imageplot(Mo);

exo4()

%% Insert your code here.
