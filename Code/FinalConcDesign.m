function [K] = FinalConcDesign

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

%% Enter the function you want to find the Fourier coefficients for:
%first for x>=0

alpha = 2/3 + 2*pi^2; beta = 2/5 + 4*pi^2 + 2*pi^4; theta = sqrt((3*(2/5 + 4*pi^2+ 2*pi^4)^2)/(2 + 6*pi^2));
NC = 5*(pi^3 + 3*pi^5)/(sqrt(6+18*pi^2)*(1+10*pi^2+5*pi^4)); NC = (1/2)/(NC);

% fx = sin(x).*(-pi<=x & x<-pi/2) + cos(3*x/2).*(-pi/2 < x & x<pi/4) + sin(x).*(pi/4 < x & x<=pi) ; %-cos(4*x).*(-pi<= x & x <=-pi/8 | 0<= x & x <= pi); 
% figure;
% plot(x,fx);
% fx = (FineX>=0).*((NC*((pi-FineX)/(sqrt(alpha)) + (pi-FineX).^3/(theta) - (beta/alpha)*(pi-FineX)/(theta)))) ... 
%     + (FineX<0).*((NC*(-pi-FineX)./(sqrt(alpha)) + (-pi-FineX).^3/(theta) - (beta/alpha)*(-pi-FineX)/(theta)));
% fx(FineX==0) = NaN;

fx_pos = -(exp(-pi)/(25*4*3*2))*(x-pi).^5;
fx_neg = (exp(-pi)/(25*4*3*2))*(-x-pi).^5;
fx = fx_pos.*(x>0).*(25*4*3*2)/(exp(-pi)*2*(pi^5)) + fx_neg.*(x<0)*(25*4*3*2)/(exp(-pi)*2*(pi^5));
fx(x==0) = 0;

figure; plot(x,fx); grid on;
title('Nonlinear Underlying Function');


fHat_1 = zeros(2*N+1,1);
for a = 1:2*N+1
fHat_1(a) = integral(@(x) (25*4*3*2)/(exp(-pi)*2*(pi^5))* -(exp(-pi)/(25*4*3*2))*(x-pi).^5 .* (exp(-1i*k(a)*x)),0,pi,'AbsTol',1e-12 );
end

fHat_2 = zeros(2*N+1,1);
for a = 1:2*N+1
fHat_2(a) = integral(@(x) (25*4*3*2)/(exp(-pi)*2*(pi^5))* (exp(-pi)/(25*4*3*2))*(-x-pi).^5 .* (exp(-1i*k(a)*x)),-pi,0,'AbsTol',1e-12 );
end

% fHat_2 = zeros(2*N+1,1);
% for a = 1:2*N+1
% fHat_2(a) = integral(@(x) (25*4*3*2)/(exp(-pi)*2*(pi^3))*(exp(-pi)/(25*4*3*2)*(-x-pi).^5).* (exp(-1i*k(a)*x)),-pi,0,'AbsTol',1e-1200 );
% end

% fx_Bernouilli = (FineX>=0).*(-1/12).*((FineX.^3)./2 + 1.5.*FineX.^2+FineX) + (-1/12).*(FineX<0).*(((FineX.^3)./2 - 1.5.*FineX.^2+FineX));
% figure; 
% plot(FineX,fx_Bernouilli); grid on; title('Bernoulli Polynomials, V_2')
% 
% fHat_1 = zeros(2*N+1,1);
% for a = 1:2*N+1
% fHat_1(a) = integral(@(FineX) (-1/12).*((FineX.^3)./2 + 1.5.*FineX.^2+FineX).* (exp(-1i*k(a)*FineX)),0,pi,'AbsTol',1e-12 );
% end
% 
% fHat_2 = zeros(2*N+1,1);
% for a = 1:2*N+1
% fHat_2(a) = integral(@(FineX) (-1/12).*((FineX.^3)./2 - 1.5.*FineX.^2+FineX).* (exp(-1i*k(a)*FineX)),-pi,0,'AbsTol',1e-12 );
% end

%fx =  (FineX>=0).*((NC*((pi-FineX)/(sqrt(alpha)) + (pi-FineX).^3/(theta) - (beta/alpha)*(pi-FineX)/(theta)))) ... 
   % + (FineX<0).*((NC*(-pi-FineX)./(sqrt(alpha)) + (-pi-FineX).^3/(theta) - (beta/alpha)*(-pi-FineX)/(theta)));

% Calculate the positive side's Fourier coefficients
% fHat_1 = zeros(2*N+1,1);
% for a = 1:2*N+1
% fHat_1(a) = integral(@(FineX) (FineX>=0).*((NC*(pi-FineX)./(sqrt(alpha))+ (pi-FineX).^3/(theta) - (beta/alpha)*(pi-FineX)/(theta))).*exp(-1i*k(a)*FineX),0,pi );
% end
% 
% % Calculate the negative side's Fourier coefficients
% fHat_2 = zeros(2*N+1,1);
% for a = 1:2*N+1
% fHat_2(a) = integral(@(FineX) (FineX<0).*((NC*(-pi-FineX)./(sqrt(alpha))+ (-pi-FineX).^3/(theta) - (beta/alpha)*(-pi-FineX)/(theta))).*exp(-1i*k(a)*FineX),-pi,0 );
% end
% 
% fHat = (fHat_1 + fHat_2)./(2*pi);
% Fourier_Approx = fHat'*exp(1i*FineX*k)';
% figure;
% plot(FineX,fx,FineX,real(Fourier_Approx));
% title('Fourier Approximation');
% 
% figure;
% plot(k,fHat); title('Fourier Coefficients');
%% Calculate the Fourier coefficients for the Gram-Schmidt orthonormalization
%  of the basis {(pi-x), (pi-x)^3}.

% 
% kern = exp(-1i*x*k);
% alpha = 2/3 + 2*pi^2; beta = 2/5 + 4*pi^2 + 2*pi^4; theta = sqrt((3*(2/5 + 4*pi^2+ 2*pi^4)^2)/(2 + 6*pi^2));
% NC = 5*(pi^3 + 3*pi^5)/(sqrt(6+18*pi^2)*(1+10*pi^2+5*pi^4)); NC = (1/2)/(NC);
% % 
% fx = (NC*((pi-FineX)/(sqrt(alpha)) + (pi-FineX).^3/(theta) - (beta/alpha)*(pi-FineX)/(theta))).*(FineX>0) + ...
%      (NC*((-pi-FineX)/(sqrt(alpha)) + (-pi-FineX).^3/(theta) - (beta/alpha)*(-pi-FineX)/(theta))).*(FineX<0);
% fx(FineX ==0) = NaN; 
 
% fx(129) = NaN;
% figure; 
% plot(x,fx); grid on;
% title('f(x),[-\pi,\pi]');
% 
% A = -(1i*pi.*k + exp(-1i*pi.*k) -1)./(k.^2); A(65) = pi^2/2; A = (A*NC)/(sqrt(alpha)*2*pi);
% B = (-1i.*pi.*k.*(-6 + pi.*k.*(pi.*k+3.*1i)) + 6.*exp(1i.*pi.*k)-6)./(k.^4); 
% B(65) = pi^4/4; B = (NC * B)/(theta * 2*pi);
% C = -(1i*pi.*k + exp(-1i*pi.*k) -1)./(k.^2); C(65) = pi^2/2; C = (C*NC*beta)/(alpha*theta*2*pi);
% 
% D = (-1i*pi.*k + exp(1i*pi.*k) -1)./(k.^2); D(65) = -(pi^2/2); D = (D*NC)/(sqrt(alpha)*2*pi);
% E = -(1i.*pi.*k.*(-6 + pi.*k.*(pi.*k-3.*1i)) + 6.*exp(1i.*pi.*k)-6)./(k.^4);
% E(65) = -(pi^4/4); E = (NC * E)/(theta * 2*pi);
% F = (-1i*pi.*k + exp(1i*pi.*k) -1)./(k.^2); F(65) = -(pi^2/2); F = (F*NC*beta)/(alpha*theta*2*pi);
% 
% fHat = A + B - C + D + E - F;
% Approx = real(fHat * kern'); 
% figure;
% plot(x,fx,x,Approx);  grid on; 
% title('Fourier Approximation, N = 64')
% legend('F(x)','Fourier Approximation')

%fx = cos(FineX).*(-pi/8 < FineX & FineX < 0) -cos(4*FineX).*((-pi <=FineX & FineX<-pi/8) | (0<FineX & FineX< pi));
% figure;
% plot(FineX,fx,'k'); grid on; title('Nonlinear Underlying Function')
% 
% %.*(-pi/8 < FineX & FineX < 0) - cos(4*FineX).*((-pi <=FineX & FineX<=-pi/8) | (0<FineX & FineX< pi)).
% 
% %Calculate the positive side's Fourier coefficients
% fHat_1 = zeros(2*N+1,1);
% for a = 1:2*N+1
% fHat_1(a) = integral(@(FineX) (NC.*((pi-FineX)/(sqrt(alpha)) + (pi-FineX).^3/(theta) - (beta/alpha)*(pi-FineX)/(theta))).*exp(-1i*k(a)*FineX),0,pi,'AbsTol',1e-12 );
% end
% 
% 
% % Calculate the negative side's Fourier coefficients
% fHat_2 = zeros(2*N+1,1);
% for a = 1:2*N+1
% fHat_2(a) = integral(@(FineX) (NC*((-pi-FineX)/(sqrt(alpha)) + (-pi-FineX).^3/(theta) - (beta/alpha)*(-pi-FineX)/(theta))).*exp(-1i*k(a)*FineX),-pi,0,'AbsTol',1e-12 );
%end

% fHat_3 = zeros(2*N+1,1);
% for a = 1:2*N+1
% fHat_3(a) = integral(@(FineX) cos(FineX).*(-pi/8 < FineX & FineX < 0) -cos(4*FineX).*((-pi <=FineX & FineX<-pi/8) | (0<FineX & FineX< pi)).*exp(-1i*k(a)*FineX),0,pi,'AbsTol',1e-12 );
% end

fHat = (fHat_1 + fHat_2)./(2*pi);
Fourier_Approx = exp(1i*x*k)*fHat;
figure; 
plot(x,fx,x,real(Fourier_Approx)); grid on; legend('g(x)','Fourier Approximation');
title('Fourier Approximation');
%% Design the new concentration factors using these Fourier coefficients. 

%fHat = 1./(2*pi*1i*k); fHat(129)=0;
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

p = 1; % Order of the factor
Concentration_poly = @(x) (p.*pi.*x.^p);
% Exponential concentration factor
alpha = 6; % Order of the factor
normal = @(x) (exp(1./(alpha*x.*(x-1))));
C = (pi)./(integral(normal,(1/N),1-(1/N)));
Concentration_exp = @(x) (C*x*exp(1./(alpha.*x.*(x-1))));

%We use the trigonometric/Gibbs concentration factor
Sipi = 1.85193705198247;					% Value of Si(pi), the sine 
                                            %integral
cfac_t = pi * sin(pi*abs(k)/N)/Sipi; 

cfe = zeros(2*N+1,1);    
cfp = zeros(2*N+1,1);
cft = zeros(2*N+1,1);
for i = 1:2*N+1
   s = abs(k(i))./(N);
   cfe(i) = Concentration_exp(s);
   cfp(i) = Concentration_poly(s);
end
% % Set cfe(n/2) and cfe(-n/2) = 0
cfe(1) = 0;
cfe(length(cfe)) = 0;

figure; plot(k,sig,k,cfe,k,cfac_t,k,cfp); grid on;
title('Concentration Factors'); legend('L1','Exponential','Trigonometric','Polynomial')

end