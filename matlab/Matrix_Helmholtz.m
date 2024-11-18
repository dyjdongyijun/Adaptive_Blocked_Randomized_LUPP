function A = Matrix_Helmholtz(NT_approx,NS_approx)

% This function gives the interaction matrix between "source" and "target" regions of different
% geometries (3-D Helmholtz)

kh     = 5.5;   % Helmholtz parameter
rad    = 3;   % Target sphere has radius "rad".
a_in   = 1;   % Inner cube is of size 2*a_in x 2*a_in x 2*a_in
%n      = ceil((NS_approx)^(1/3));    % Number of Chebyshev nodes along a side of inner box.
n = 20;

%%% Set up the source geometry (cube)
[x,w]      = LOCAL_clencurt(n,-a_in,a_in);
x          = x(end:(-1):1)';
w          = w';
[XX1,XX2,XX3]  = meshgrid(x);
xxs        = [reshape(XX1,1,n*n*n);...
              reshape(XX2,1,n*n*n);...
              reshape(XX3,1,n*n*n)];
tmp = zeros(length(w),length(w),length(w));
for ii = 1:length(w)
    for jj = 1:length(w)
        for kk = 1:length(w)
            tmp(ii,jj,kk) = w(ii)*w(jj)*w(kk);
        end
    end
end
wws        = reshape(tmp,1,n*n*n);
NS         = size(xxs,2);

%%% Set up the target geometry (sphere)
[xxt(1,:),xxt(2,:),xxt(3,:),NT] = mySphere(NT_approx);
xxt = rad.*xxt;

%%% Construct NT x NS interaction matrix
Afun = @(i,j)Afun_(i,j,xxt,xxs,kh,wws);
A = Afun(1:NT,1:NS);

end

function [x,w] = LOCAL_clencurt(N1,a,b)

% Adapted from "fclencurt.m" by Greg von Winckel.
N=N1-1; 
bma=b-a;
c=zeros(N1,2);
c(1:2:N1,1)=(2./[1 1-(2:2:N).^2 ])'; 
c(2,2)=1;
f=real(ifft([c(1:N1,:);c(N:-1:2,:)]));
w=bma*([f(1,1); 2*f(2:N,1); f(N1,1)])/2;
x=0.5*((b+a)+N*bma*f(1:N1,2));

return

end

function [X,Y,Z,N_new] = mySphere(N)
%% Generate Node xyz positions
% Used 2004 paper by Markus Deserno, Max-Planck-Institut:
% "How to generate equidistributed points on the surface of a sphere"
% Enforces constant intervales d_theta ~ d_phi
% Assumes unit radius
% Does not replace MATLAB "sphere" function
% Create Sphere 3D Geometry Centered at (x,y,z) = (0,0,0)
%
% N: target number of nodes
% N_new: final number of nodes
% X,Y,Z: column vectors of length N_new containing node coordinates
disp("In mySphere");
r_unit = 1;
Area = 4*pi*r_unit^2/N;
Distance = sqrt(Area);
M_theta = round(pi/Distance);
d_theta = pi/M_theta;
d_phi = Area/d_theta;
N_new = 0;
X = nan(N_new,1);
Y = nan(N_new,1);
Z = nan(N_new,1);
disp(M_theta);
for m = 0:M_theta-1
    
    Theta = pi*(m+0.5)/M_theta;
    M_phi = round(2*pi*sin(Theta)/d_phi); % not exact
    
    for n = 0:M_phi-1        
        Phi = 2*pi*n/M_phi;    
        
        N_new = N_new + 1;
        
        X(N_new) = sin(Theta)*cos(Phi);
        Y(N_new) = sin(Theta)*sin(Phi);
        Z(N_new) = cos(Theta);
        
    end
end
end

% kernel function
function K = Kfun(x,y,k,w)
  dx = x(1,:)' - y(1,:);
  dy = x(2,:)' - y(2,:);
  dz = x(3,:)' - y(3,:);
  dr = sqrt(dx.^2 + dy.^2 + dz.^2);
  NT = size(x,2);
  K = (1/(4*pi)*exp(1i*k*dr)./dr).*(ones(NT,1) * w);
end

% matrix entries
function A = Afun_(i,j,rx,cx,k,w)
  A = Kfun(rx(:,i),cx(:,j),k,w(:,j));
end

