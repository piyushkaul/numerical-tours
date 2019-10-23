
warning off
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/denoisingwav_3_advanced')
warning on

n = 256;
name = 'hibiscus';
M0 = load_image(name,n);
M0 = rescale( sum(M0,3) );

sigma = .08;

M = M0 + sigma*randn(size(M0));

Jmin = 3;
MW = perform_wavelet_transf(M,Jmin,+1);

T = 1; % threshold value
v = linspace(-5,5,1024);
clf;
hold('on');
plot(v, perform_thresholding(v,T,'hard'), 'b--');
plot(v, perform_thresholding(v,T,'soft'), 'r--');
plot(v, perform_thresholding(v,[T 2*T],'semisoft'), 'g');
plot(v, perform_thresholding(v,[T 4*T],'semisoft'), 'g:');
legend('hard', 'soft', 'semisoft, \mu=2', 'semisoft, \mu=4');
hold('off');

exo1()

%% Insert your code here.

err_mu = compute_min(err, 2);
clf;
plot(mulist, err_mu, '.-');
axis('tight');
set_label('\mu', 'SNR');

T = 1; % threshold value
v = linspace(-4,4,1024);

v_hard = v .* (abs(v)>T);

v_soft = v .* max( 1-T./abs(v), 0 );

v_stein = v .* max( 1-(T^2)./(v.^2), 0 );

clf;
hold('on');
plot(v, v_hard, 'b');
plot(v, v_soft, 'r');
plot(v, v_stein, 'k--');
hold('off');
legend('Hard', 'Soft', 'Stein');

exo2()

%% Insert your code here.
