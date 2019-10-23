
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('toolbox_graph')
addpath('solutions/graphics_4_multiplicative_cascade')

T = 10;

rmin = 0.02;

Delta_t = rmin/2;
n = T/Delta_t+1;

lambda = (1/rmin-1)*(T+1); % to ensure scale invariance, scales are distributed as 1/r^2.
N = round(lambda);

ti = -1/2+rand(1,N) * (T+1);
ti_1 = -1/2+rand(1,round(T+1))*(T+1);
ti = [ti_1 ti];   % for exact scale invariance

umax = 1/rmin-1;
ui = [zeros(1,length(ti_1)) rand(1,N) * umax ]; % ui = 1/ri-1
ri = (1+ui).^(-1);

figure(1)
clf
plot(ti, ri, '.')
axis([-1/2,T+1/2,0,1]);

sigma2 = 0.2;
mu = -sigma2/2;

if -(exp(2*(sigma2+mu))-1)<-1
   disp('Be careful ! This cascade will degenerate as rmin -> 0 !')
end

N = length(ti);   % the number of multipliers = number of time-scale points.
Wi = exp( randn(N,1)*sqrt(sigma2)+mu );

t = linspace(0,T,n);

H1 = 1 - exp(mu+sigma2/2);
f = ones(n,1) * exp(H1) / rmin^H1;

i = 1;
I = find(abs(t-ti(i))<=ri(i)/2); % t belongs to a disk centered on ti(i)

f(I) = f(I) * Wi(i);

clf
plot(t,f)
axis([0 T 0 1.1*max(f)])

exo1()

%% Insert your code here.

figure(1)
clf
plot(t,f)
axis([0 T 0 1.1*max(f)])

exo2()

%% Insert your code here.

exo3()

%% Insert your code here.

n = 128;

rmin = 1/n;

Xmax = 1;
Ymax = 1;

lambda = 2/pi*(1/rmin^2-1)*(Xmax+1)*(Ymax+1); % density g(r)dr=4/pi/r^3 dr
N = round(lambda); % should be a Poisson r.v. with expectation lambda.

umax = (1/rmin^2-1);          % u will be a uniform variable in [0 1/rmin^2-1]
ui = rand(1,N) * umax;
ri = (1/rmin^2-ui).^(-1/2);

xi = -1/2 + rand(1,N) * (Xmax+1);
yi = -1/2 + rand(1,N) * (Ymax+1);

clf
h = plot3(xi, yi, ri, '.');
axis([-1/2,Xmax+1/2,-1/2,Ymax+1/2,0,1]);

sigma2 = 0.08;
mu = -sigma2/2;

Wi = exp( randn(N,1)*sqrt(sigma2)+mu );

x = linspace(0,Xmax,n);
y = linspace(0,Ymax,n);
[X,Y]= meshgrid(x,y);

H1 = 1 - exp(mu+sigma2/2);
f = ones(n)/rmin^H1;

i = 1;
I = find( (X-xi(i)).^2+(Y-yi(i)).^2 <=ri(i)^2/4 );

f(I) = f(I) * Wi(i);

exo4()

%% Insert your code here.

clf
imageplot(f)

x = [0:n/2 -n/2+1:-1];
[U,V] = meshgrid(x,x); 
S = sqrt(U.^2+V.^2); 
S(1,1) = 1;

alpha = .5;

F = real( ifft2( fft2(f)./S.^alpha ) );

exo5()

%% Insert your code here.

exo6()

%% Insert your code here.
