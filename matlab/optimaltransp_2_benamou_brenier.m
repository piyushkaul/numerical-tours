
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/optimaltransp_2_benamou_brenier')

mynorm = @(a)norm(a(:));
sum3 = @(a)sum(a(:));

n = 20;
p = 20;

[Y,X] = meshgrid(linspace(0,1,n), linspace(0,1,n));
gaussian = @(a,b,sigma)exp( -((X-a).^2+(Y-b).^2)/(2*sigma^2) );
normalize = @(u)u/sum(u(:));

sigma = .1;
rho = .05; % minimum density value
f0 = normalize( rho + gaussian(.2,.3,sigma) );
f1 = normalize( rho + gaussian(.6,.7,sigma*.7) + .6*gaussian(.7,.4,sigma*.7) );

clf;
imageplot({f0 f1});

bound = 'per';
bound = 'neum';

if strcmp(bound, 'per')
    dx = @(u)u([2:end 1],:,:)-u;
    dy = @(u)u(:,[2:end 1],:)-u;
else
    dx = @(u)u([2:end end],:,:)-u;
    dy = @(u)u(:,[2:end end],:)-u;
end

if strcmp(bound, 'per')
    dxS = @(u)-u+u([end 1:end-1],:,:);
    dyS = @(u)-u+u(:,[end 1:end-1],:);
else
    dxS = @(u)[-u(1,:,:); u(1:end-2,:,:)-u(2:end-1,:,:); u(end-1,:,:)];
    dyS = @(u)[-u(:,1,:), u(:,1:end-2,:)-u(:,2:end-1,:), u(:,end-1,:)];    
end

fprintf('Should be 0: %.2e\n', certify_adjoint(dx,dxS,[n n p]));
fprintf('Should be 0: %.2e\n', certify_adjoint(dy,dyS,[n n p]));

grad = @(f)cat(4, dx(f), dy(f));
div  = @(u)-dxS(u(:,:,:,1)) - dyS(u(:,:,:,2));

fprintf('Should be 0: %.2e\n', certify_adjoint(grad,@(v)-div(v),[n n p]));

dt  = @(f)cat(3, f(:,:,2:end)-f(:,:,1:end-1), zeros(size(f,1),size(f,2)) );
dtS = @(u)cat(3, -u(:,:,1), u(:,:,1:end-2)-u(:,:,2:end-1), u(:,:,end-1));

A = @(w)cat( 3, div(w(:,:,:,1:2))+dt(w(:,:,:,3)), w(:,:,1,3), w(:,:,end,3) );

U = @(r0,r1)cat(3, r0, zeros(n,n,p-2), r1);
AS = @(s)cat(4, -grad(s(:,:,1:p)), dtS(s(:,:,1:p)) + U(s(:,:,end-1),s(:,:,end)) );

fprintf('Should be 0: %.2e\n', certify_adjoint(A,AS,[n n p 3]));

r0 = cat(3, zeros(n,n,p), f0, f1);

J = @(w)sum3(  sum(w(:,:,:,1:2).^2,4) ./ w(:,:,:,3)   );

PolyCoef = @(m0,f0,lambda)[ones(length(f0),1), 4*lambda-f0, 4*lambda^2-4*f0, -lambda*sum(m0.^2,2) - 4*lambda^2*f0];

extract = @(A)A(:,1);
CubicReal = @(P)real( extract(poly_root(P')') );

Proxj0 = @(m0,f, lambda)cat(2, m0 ./ repmat( 1+2*lambda./f, [1 2]), f );
Proxj  = @(m0,f0,lambda)Proxj0( m0, CubicReal(PolyCoef(m0,f0,lambda)), lambda );

ProxJ = @(w,lambda)reshape( Proxj( ...
                   reshape(w(:,:,:,1:2), [n*n*p 2]), ...
                   reshape(w(:,:,:,3  ), [n*n*p 1]), lambda ), [n n p 3] );

opts.epsilon = 1e-9; 
opts.niter_max = 150;

flat = @(x)x(:);
resh = @(x)reshape(x, [n n p+2]);
mycg = @(B,y)resh( perform_cg(@(r)flat(B(resh(r))),y(:),opts) );

pA = @(r)mycg(@(s)A(AS(s)),r);

ProxG = @(w,lambda)w + AS( pA(r0-A(w)) );

w = randn(n,n,p,3);
err = @(w)mynorm(A(w)-r0)/mynorm(r0);
fprintf('Error before projection: %.2e\n', err(w));
fprintf('Error before projection: %.2e\n', err(ProxG(w)));

mu = 1;
gamma = 1;

rProxJ = @(w,tau)2*ProxJ(w,tau)-w;
rProxG = @(w,tau)2*ProxG(w,tau)-w;

niter = 200;

t = repmat( reshape(linspace(0,1,p), [1 1 p]), [n n 1]);
f = (1-t) .* repmat(f0, [1 1 p]) + t .* repmat(f1, [1 1 p]);
m = zeros(n,n,p,2);
w0 = cat(4, m,f);

sel = round(linspace(1,p,6));
clf;
imageplot( mat2cell(w0(:,:,sel,3), n, n, ones(6,1)) , '', 2,3);

w = ProxG(w0,gamma);
mynorm(A(w0)-r0)/mynorm(r0)
mynorm(A(w)-r0)/mynorm(r0)

exo1()

%% Insert your code here.

sel = round(linspace(1,p,6));
clf;
imageplot( mat2cell(w(:,:,sel,3), n, n, ones(6,1)) , '', 2,3);
