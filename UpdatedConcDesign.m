
%% Set up Fourier modes and physical mesh.
% Clears all variables in the workspace and closes all open figures
clear all; close all; clc

% Alter figure label and text font sizes
set(0,'defaultaxesfontsize',15);
set(0,'defaulttextfontsize',15); 

N = 64;           % No. of modes to plot is 2N+1
k = (-N:N);       % The Fourier Modes
    
% The physical grid that is used for generating the concentration factors.
nPts = 256;
h = 2*pi/nPts;
x = -pi + h*(0:nPts-1).';


%% Define the function we are going to use to derive the new
%  new concentration factors.

% The template function comes from a Gram-Schmidt orhtonormalization of the
% basis {(x-pi),(x-pi)^3}, normalized to have a jump of 1 at x=0. 


kern = exp(-1i*x*k);
kern_2 = exp(1i*x*k); 
alpha = 2/3 + 2*pi^2; beta = 2/5 + 4*pi^2 + 2*pi^4; theta = sqrt((3*(2/5 + 4*pi^2+ 2*pi^4)^2)/(2 + 6*pi^2));
NC = 5*(pi^3 + 3*pi^5)/(sqrt(6+18*pi^2)*(1+10*pi^2+5*pi^4)); NC = (1/2)/(NC);

fx = (NC*((pi-x)/(sqrt(alpha)) + (pi-x).^3/(theta) - (beta/alpha)*(pi-x)/(theta))).*(x>0) + ...
    (NC*((-pi-x)/(sqrt(alpha)) + (-pi-x).^3/(theta) - (beta/alpha)*(-pi-x)/(theta))).*(x<0);
fx(129) = NaN;
figure; 
plot(x,fx); grid on;
title('f(x),[-\pi,\pi]')

A = -(1i*pi.*k + exp(-1i*pi.*k) -1)./(k.^2); 
A(65) = pi^2/2; A = (A*NC)/(sqrt(alpha)*2*pi);
B = (-1i.*pi.*k.*(-6 + pi.*k.*(pi.*k+3.*1i)) + 6.*exp(1i.*pi.*k)-6)./(k.^4);
B(65) = pi^4/4; B = (NC * B)/(theta * 2*pi);
C = -(1i*pi.*k + exp(-1i*pi.*k) -1)./(k.^2); C(65) = pi^2/2; C = (C*NC*beta)/(alpha*theta*2*pi);

D = (-1i*pi.*k + exp(1i*pi.*k) -1)./(k.^2); D(65) = -(pi^2/2); D = (D*NC)/(sqrt(alpha)*2*pi);
E = -(1i.*pi.*k.*(-6 + pi.*k.*(pi.*k-3.*1i)) + 6.*exp(1i.*pi.*k)-6)./(k.^4);
E(65) = -(pi^4/4); E = (NC * E)/(theta * 2*pi);
F = (-1i*pi.*k + exp(1i*pi.*k) -1)./(k.^2); F(65) = -(pi^2/2); F = (F*NC*beta)/(alpha*theta*2*pi);

fHat = A + B - C + D + E - F;
Approx = real(fHat * kern'); 
figure;
plot(x,fx,x,Approx);  grid on; 

%% This code uses cvx to derive new concentration factors.

cvx_begin

variable sig(2*N+1) 

%variable Final_Kernel

K = zeros(length(x),2*N +1);

for j = 1:length(x)
    for a = 1:2*N +1 
         K(j,a) = 1i*fHat(a)*sign(k(a))*exp(1i*k(a) * x(j));
    end
end

minimize norm(K*sig,1)
subject to 
    K(x==0,:)*sig == 1;
cvx_end

figure; plot(k,sig); grid on; title('Concentration Factor');
figure; plot(x,real(K*sig)); grid on; title('Jump Approximation');

