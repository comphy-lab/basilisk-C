%% Study of ginsburg-landau eq in a strongly convective regime
%
% dF/dt = -U dF/dx + sigma(x) F + kappa d^2F/dx^2
%
%


xmin = -10;xmax =50;N = 500;
discretization = 'fd'; 
% try either 'fds' (finite differences), 'fd' (fd with sparse matrix storage) or  'cheb' (chebyshev)
[dx,dxx,wx,x]=dif1D(discretization,xmin,xmax-xmin,N);
I = diag(ones(N,1));Z = zeros(N,N);
symbols = {'r','b','g','c','m','y'};
symbolsS = {'*','+','^','v','o','p'};
if ~exist('isplot')
    isplot = 1; isplotS = 1;
else
    isplot = isplot+1;
    if isplot>length(symbols)
        isplot = 1; isplotS = isplotS+1;
    end
end

%% model A
%U = 2;
%sigma0=1;
%sigma1=2;
%kappa=1+1i;
%sigmax = (sigma0+sigma1./x);

%% model B
%U = 1;
%kappa = 1-i;
%sigma0 = 1;sigma2 = 1/16;
%sigmax = sigma0-sigma2*x.^2/2;

%% model C
U = 2
kappa = 1;
sigma0 = 1;sigma2 = 1/50;sigmainf = 0.5;

sigmaabsmin = real(U^2/(4*kappa))
sigmax = sigma0-sigma2*x.^2./(2+sigma2/(sigma0-sigmainf)*(x>0).*x.^2);

%% theoretical solution (WKB approx)

stheo = (sigma0-U^2/(4*kappa))*[1,1,1,1]-sqrt(2*sigma2*kappa)*[1/2,3/2,5/2,7/2];
l= (-5:.01:5);
scontinuous = sigmainf-1i*U*l-kappa*l.^2;
figure(4);plot(imag(stheo),real(stheo),'*k',imag(scontinuous),real(scontinuous),':k');hold on;

%% plot absolute growth rate
figure(1); plot(x,sigmax,'r',x,real(U^2/(4*kappa))*ones(length(x)),'b');

%% complex mapping definition
Xc = 5;Lc = 5;gammac = -1;
g = (x>Xc).*tanh((x-Xc).^2/Lc.^2);

Gx = x.*(x<Xc)+(x>=Xc).*(Xc+(1+1i*gammac*g).*(x-Xc));
Hx = 1*(x<Xc)+(x>Xc).*(0.1e1 ./ (2*1i * gammac * ((x - Xc) .^ 2) / (Lc ^ 2) .* (0.1e1 - tanh((x - Xc) .^ 2 / Lc ^ 2) .^ 2)...
    + 1 + 1i * gammac * tanh((x - Xc) .^ 2 / Lc ^ 2)));
dxM = Hx.*dx;
dxxM = dxM^2;

figure(10); plot(x,real(Gx),x,imag(Gx),x,real(Hx),x,imag(Hx));

%% build matrix
%A = -U*dx+sigmax.*I+kappa*dxx;
A = -U*dxM+sigmax.*I+kappa*dxxM;
B = I;

A(1,:)=-1e-10*I(1,:);A(N,:)=-1e-10*I(N,:);B(1,:)=Z(1,:);B(N,:)=Z(N,:);

[UU,S]=eig(A,B);

s = diag(S);
%[t,o]=sort(abs(s)); 
[t,o]=sort(-real(s)); 
s=s(o); UU=UU(:,o);


[~,ind] = min(abs(s-stheo(1)));[~,ind2] = min(abs(s-stheo(2)));
figure(3); plot(x,UU(:,ind),[symbols{isplot},symbolsS{isplotS},'-']);hold on;
hold on; plot(x,UU(:,ind2),[symbols{isplot},symbolsS{isplotS},'--']);


figure(4);hold on; plot(imag(s),real(s),[symbols{isplot},symbolsS{isplotS}]);
xlim([-10 20]);ylim([-10 2]);

%%
scontinuousT = sigmainf-1i*U*(l./(1+1i*gammac))-kappa*(l.^2./(1+1i*gammac)^2);
figure(4); plot(imag(scontinuousT),real(scontinuousT),[symbols{isplot},':']);


%% compute the resolvant
omegatab = 0:.01:.5;
for jj = 1:length(omegatab)
    omega = omegatab(jj);
    Res(jj) = norm(inv(A-1i*omega*B));
end
figure(40); hold on; plot(omegatab,Res,[symbols{isplot},'-']);

%%    
omega = 0.36;
[UUU,S,V] = svd(inv(A-1i*omega*B));
figure(33); plot(x,UUU(:,1),[symbols{isplot},'-']);
