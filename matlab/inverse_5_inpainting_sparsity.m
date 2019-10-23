
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/inverse_5_inpainting_sparsity')

n = 128;
name = 'lena';
f0 = load_image(name);
f0 = rescale(crop(f0,n));

clf;
imageplot(f0, 'Image f_0');

rho = .7;

Omega = zeros(n,n);
sel = randperm(n^2); 
Omega(sel(1:round(rho*n^2))) = 1;

if using_matlab()
    Phi = @(f,Omega)f.*(1-Omega);
end

y = Phi(f0,Omega);

clf;
imageplot(y, 'Observations y');

if using_matlab()
    SoftThresh = @(x,T)x.*max( 0, 1-T./max(abs(x),1e-10) );
end

clf;
T = linspace(-1,1,1000);
plot( T, SoftThresh(T,.5) );
axis('equal');

Jmax = log2(n)-1;
Jmin = Jmax-3;

options.ti = 0; % use orthogonality.
if using_matlab()
    Psi = @(a)perform_wavelet_transf(a, Jmin, -1,options);
    PsiS = @(f)perform_wavelet_transf(f, Jmin, +1,options);
end

if using_matlab()
    SoftThreshPsi = @(f,T)Psi(SoftThresh(PsiS(f),T));
end

clf;
imageplot( clamp(SoftThreshPsi(f0,.1)) );

lambda = .03;

if using_matlab()
    ProjC = @(f,Omega)Omega.*f + (1-Omega).*y;
end

fSpars = y;

fSpars = ProjC(fSpars,Omega);

fSpars = SoftThreshPsi( fSpars, lambda );

exo1()

%% Insert your code here.

clf;
imageplot(clamp(fSpars));

exo2()

%% Insert your code here.

J = Jmax-Jmin+1;
u = [4^(-J) 4.^(-floor(J+2/3:-1/3:1)) ];
U = repmat( reshape(u,[1 1 length(u)]), [n n 1] );

lambda = .01;

options.ti = 1; % use translation invariance
if using_matlab()
    Xi = @(a)perform_wavelet_transf(a, Jmin, -1,options);
    PsiS = @(f)perform_wavelet_transf(f, Jmin, +1,options);
    Psi = @(a)Xi(a./U);
end

tau = 1.9*min(u);

a = U.*PsiS(fSpars);

fTI = Psi(a);    
a = a + tau*PsiS( Phi( y-Phi(fTI,Omega),Omega ) );

a = SoftThresh( a, lambda*tau );

exo3()

%% Insert your code here.

fTI = Psi(a);

clf;
imageplot(clamp(fTI));

exo4()

%% Insert your code here.

if using_matlab()
    HardThresh = @(x,t)x.*(abs(x)>t);
end

t = linspace(-1,1,1000);
plot( t, HardThresh(t,.5) );
axis('equal');

niter = 500;

lambda_list = linspace(1,0,niter);

fHard = y;

fHard = ProjC(fHard,Omega);

fHard = Xi( HardThresh( PsiS(fHard), tau*lambda_list(1) ) );

exo5()

%% Insert your code here.
