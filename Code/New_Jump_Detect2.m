%% Receiver Operating Characteristic (ROC) curve
% This script plots the ROC for the Gibbs concentration factor 
% and a 3-point detector. With slight adjustments the code can be made to
% run with a 5-point detection vector.

% change the resolution of the fineX grid, make them the same as the FC.
% Make sure they compare to the right threshold.

close all; clear all; clc
%% Initialization 
% Define the Fourier modes
nModes = 64;
fourModes = (-nModes:nModes).';
k = (-nModes:nModes)';

% Generate a (periodic) physical space grid
nGridPts = 256;
xl = -pi; xr = pi; h = (xr-xl)/nGridPts;
x = xl + h*(0:nGridPts-1).';
%% Here we are going to design the concentration factors that will be used 
% to run the minmod tests, the third threshold test.
% Polynomial concentration factor
p = 1; % Order of the factor
Concentration_poly = @(x) (p.*pi.*x.^p);
% Exponential concentration factor
alpha = 6; % Order of the factor
normal = @(x) (exp(1./(alpha*x.*(x-1))));
C = (pi)./(integral(normal,(1/nModes),1-(1/nModes)));
Concentration_exp = @(x) (C*x*exp(1./(alpha.*x.*(x-1))));

%We use the trigonometric/Gibbs concentration factor
Sipi = 1.85193705198247;					% Value of Si(pi), the sine 
                                            %integral
cfac_t = pi * sin(pi*abs(fourModes)/nModes)/Sipi; 

cfe = zeros(2*nModes+1,1);    
cfp = zeros(2*nModes+1,1);
cft = zeros(2*nModes+1,1);
for i = 1:2*nModes+1
   s = abs(k(i))./(nModes);
   cfe(i) = Concentration_exp(s);
   cfp(i) = Concentration_poly(s);
end
% % Set cfe(n/2) and cfe(-n/2) = 0
cfe(1) = 0;
cfe(length(cfe)) = 0;

figure(1);
plot(k,cfe,'r',k,cfp,'g',k,cfac_t,'b');
title('Concentration Factors, N = 128'); grid on;
legend('Exponential', 'Polynomial', 'Trigonometric');
%% Concentration factor definitions
%We use the trigonometric/Gibbs concentration factor
Sipi = 1.85193705198247;					% Value of Si(pi), the sine 
                                            %integral
cfac_t = pi * sin(pi*abs(fourModes)/nModes)/Sipi;   
sig = UpdatedConcDesign;
 
cfac_1 = cfac_t;
cfac_2 = sig;
cfac_3 = sig;
cfac_4 = cfac_t;
cfac_5 = sig;
figure;
plot(k,sig); title('Concentration Factor'); grid on; 

%% Here we calculate the Fourier coefficients for a function, using a fine x-mesh.
nPts = 256;
h = 2*pi/nPts;
FineX = -pi + h*(0:nPts-1).';
N = nModes;

alpha = 2/3 + 2*pi^2; beta = 2/5 + 4*pi^2 + 2*pi^4; theta = sqrt((3*(2/5 + 4*pi^2+ 2*pi^4)^2)/(2 + 6*pi^2));
NC = 5*(pi^3 + 3*pi^5)/(sqrt(6+18*pi^2)*(1+10*pi^2+5*pi^4)); NC = (1/2)/(NC);

% figure;
% plot(x,fx);
fx = (x>=0).*((NC*((pi-x)/(sqrt(alpha)) + (pi-x).^3/(theta) - (beta/alpha)*(pi-x)/(theta)))) ... 
    + (x<0).*((NC*(-pi-x)./(sqrt(alpha)) + (-pi-x).^3/(theta) - (beta/alpha)*(-pi-x)/(theta)));


%fx =  (FineX>=0).*((NC*((pi-FineX)/(sqrt(alpha)) + (pi-FineX).^3/(theta) - (beta/alpha)*(pi-FineX)/(theta)))) ... 
   % + (FineX<0).*((NC*(-pi-FineX)./(sqrt(alpha)) + (-pi-FineX).^3/(theta) - (beta/alpha)*(-pi-FineX)/(theta)));

%fx = 3/2.*(-.75*pi <= x & x< -pi/2) + (7/4 -.5*x + sin(x-.25)).*(-.25*pi <x & x< pi/8) + (11*x/4-5).*(3*pi/8 <=x & x<3*pi/4);
%fx  = cos(3*x/2).*(-pi/8 < x & x<0) + cos(11*x).*(-pi <= x & x<= -pi/8 | 0<= x <=pi);
%fx = cos(3*x/2).*(abs(x)< pi/2) + cos(7*x/2).*(abs(x)>=pi/2);
% %Calculate the positive side's Fourier coefficients
% fHat_1 = zeros(2*N+1,1);
% for a = 1:2*N+1
% fHat_1(a) = integral(@(FineX) ((NC*((pi-x)/(sqrt(alpha)) + (pi-x).^3/(theta) - (beta/alpha)*(pi-x)/(theta)))).*exp(-1i*k(a)*FineX),0,pi,'AbsTol',1e-12 );
% end
% 
% 
% % Calculate the negative side's Fourier coefficients
% fHat_2 = zeros(2*N+1,1);
% for a = 1:2*N+1
% fHat_2(a) = integral(@(FineX)  ((NC*((pi-x)/(sqrt(alpha)) + (-pi-x).^3/(theta) - (beta/alpha)*(-pi-x)/(theta)))).*exp(-1i*k(a)*FineX),-pi,0,'AbsTol',1e-12 );
% end

figure;
plot(FineX,fx); grid on;
title('Nonlinear Underyling Function')
%% Use this if you want to use the function from Peterson, Gelb, and Eubank:

%Calculate the positive side's Fourier coefficients
% fHat_1 = zeros(2*N+1,1);
% for a = 1:2*N+1
% fHat_1(a) = integral(@(FineX) cos(7*FineX).*exp(-1i*k(a)*FineX),-pi/8,0,'AbsTol',1e-12 );
% end
% 
% 
% % Calculate the negative side's Fourier coefficients
% fHat_2 = zeros(2*N+1,1);
% for a = 1:2*N+1
% fHat_2(a) = integral(@(FineX)  -cos(11*FineX).*exp(-1i*k(a)*FineX),-pi,-pi/8,'AbsTol',1e-12 );
% end
% 
% fHat_3 = zeros(2*N+1,1);
% for a = 1:2*N+1
% fHat_3(a) = integral(@(FineX) -cos(11*FineX).*exp(-1i*k(a)*FineX),0,pi,'AbsTol',1e-12 );
% end
% 
% fx = cos(FineX).*(-pi/8 < FineX & FineX < 0) -cos(4*FineX).*((-pi <=FineX & FineX<-pi/8) | (0<FineX & FineX< pi));
% figure;
%  plot(x,fx);

%% Use this is you want to use the function from Aditya's thesis:
% fHat_1 =zeros(2*N+1,1);
% for a =1:2*N+1
%    fHat_1(a) = integral(@(FineX) 1.5.*exp(-1i.*k(a).*FineX),-.75*pi,-pi/2); 
% end
% 
% fHat_2 =zeros(2*N+1,1);
% for a =1:2*N+1
%    fHat_2(a) = integral(@(FineX) (7./4 - .5*FineX + sin(FineX-.25)).*exp(-1i.*k(a).*FineX),-.25*pi,pi/8); 
% end
% 
% fHat_3 =zeros(2*N+1,1);
% for a =1:2*N+1
%    fHat_3(a) = integral(@(FineX) (11.*FineX./4 - 5).*exp(-1i*k(a)*FineX),3*pi/8,3*pi/4); 
% end
% 
% fx = (1.5 .*(x>-3*pi/4 & x<-pi/2)) + (7./4 - .5*FineX + sin(FineX-.25)).*(x>-pi/4 & x<pi/8) + (11.*FineX./4 - 5).*(x>3*pi/8 & x<3*pi/4);
% fHat_1 = zeros(2*N+1,1);
% for a = 1:2*N+1
%     fHat_1(a) = integral(@(FineX) cos(3*FineX/2).*exp(-1i*k(a)*FineX),-pi/2,pi/2,'AbsTol',1e-12);
% end
% 
% fHat_2 = zeros(2*N+1,1);
% for a = 1:2*N+1
%     fHat_2(a) = integral(@(FineX) cos(7.*FineX./2).*exp(-1i*k(a)*FineX),pi/2,pi,'AbsTol',1e-12);
% end
% 
% fHat_3 = zeros(2*N+1,1);
% for a = 1:2*N+1
%     fHat_3(a) = integral(@(FineX) cos(7.*FineX./2).*exp(-1i*k(a)*FineX),-pi,-pi/2,'AbsTol',1e-12);
% end


fHat = (fHat_1 + fHat_2+ fHat_3)./(2*pi);
Func_Approx = fHat' * (exp(1i*FineX*k'))';
figure;
plot(FineX,Func_Approx);

figure;
plot(k,real(fHat)); title('Fourier Coefficients');

%% Function definition
% The template function is a periodic ramp function with unit jump at x = 0
%fx = sin(2*x).*(x<0) + exp(-3*x).*(x>=0);
% kern = exp(-1i*x*k');
% kern_2 = exp(1i*x*k'); 
% alpha = 2/3 + 2*pi^2; beta = 2/5 + 4*pi^2 + 2*pi^4; theta = sqrt((3*(2/5 + 4*pi^2+ 2*pi^4)^2)/(2 + 6*pi^2));
% NC = 5*(pi^3 + 3*pi^5)/(sqrt(6+18*pi^2)*(1+10*pi^2+5*pi^4)); NC = (1/2)/(NC);
% 
% fx = (NC*((pi-x)/(sqrt(alpha)) + (pi-x).^3/(theta) - (beta/alpha)*(pi-x)/(theta))).*(x>0) + ...
%     (NC*((-pi-x)/(sqrt(alpha)) + (-pi-x).^3/(theta) - (beta/alpha)*(-pi-x)/(theta))).*(x<0);
% fx(129) = NaN;
% figure; 
% plot(x,fx); grid on;
% title('f(x),[-\pi,\pi]')
% 
% A = -(1i*pi.*k + exp(-1i*pi.*k) -1)./(k.^2); 
% A(65) = pi^2/2; A = (A*NC)/(sqrt(alpha)*2*pi);
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
% Approx = real(fHat' * kern_2'); 
% figure;
% plot(x,fx,x,Approx);  grid on; 
% 
% 
fourMat = exp(1i*FineX*fourModes.' );


% Template or jump response
jmpFncCfs =  1i*fHat' .* sign(fourModes)' .* cfe';
template =  fourMat*jmpFncCfs' ;
figure; 
plot(FineX,real(template)); grid on; 
%% Noise characteristics
% We consider zero mean, additive white Gaussian nose of variance rho^2
rho2 = 1/nModes^2;
n = sqrt(rho2)*randn(2*nModes+1,1) + 1i*sqrt(rho2)*randn(2*nModes+1,1);
fHat_alt  = fHat + n;

%delta
delta = 1/(2*nModes + 1);
%define a probability for false alarms
pfa = .01;

jump_detect = zeros(length(FineX),1);
test_stat_adv = zeros(length(FineX),1);
gamma_dot = NaN(length(FineX),1);
test_stat = NaN(length(FineX),1);
cov = zeros(5,5,length(FineX));
for a = 1:length(FineX)

% First calculate the nonshifted jump responses:  
    %Define stencil
    stencil(1) = FineX(a) - 2*delta;
    stencil(2) = FineX(a) - delta;
    stencil(3) = FineX(a);
    stencil(4) = FineX(a) + delta;
    stencil(5) = FineX(a) + 2*delta;
    %Define M
     M = zeros(5,1);
     Exponent_factor_1 = exp( 1i*(stencil(1)) * fourModes.');
     Exponent_factor_2 = exp( 1i*(stencil(2)) * fourModes.');
     Exponent_factor_3 = exp( 1i*(stencil(3)) * fourModes.');
     Exponent_factor_4 = exp( 1i*(stencil(4)) * fourModes.');
     Exponent_factor_5 = exp( 1i*(stencil(5)) * fourModes.');
     jmpFncCfs_1 = 1i * fHat' .* sign(fourModes)' .* cfac_1';
     jmpFncCfs_2 = 1i * fHat' .* sign(fourModes)' .* cfac_2';
     jmpFncCfs_3 = 1i * fHat' .* sign(fourModes)' .* cfac_3';
     jmpFncCfs_4 = 1i * fHat' .* sign(fourModes)' .* cfac_4';
     jmpFncCfs_5 = 1i * fHat' .* sign(fourModes)' .* cfac_5';
        M_1 = real(Exponent_factor_1 * jmpFncCfs_1');
        M_2 = real(Exponent_factor_2 * jmpFncCfs_2');
        M_3 = real(Exponent_factor_3 * jmpFncCfs_3');
        M_4 = real(Exponent_factor_4 * jmpFncCfs_4');
        M_5 = real(Exponent_factor_5 * jmpFncCfs_5');
        M(1) = M_1; M(2) = M_2; M(3) = M_3; 
        M(4) = M_4; M(5) = M_5;
        M = M;
    %Define covariance matrix

     cov(1,1,a) = sum((rho2) * real(cfac_1.^2));
     cov(2,2,a) = sum((rho2) * real(cfac_2.^2));
     cov(3,3,a) = sum((rho2) * real(cfac_3.^2));
     cov(4,4,a) = sum((rho2) * real(cfac_4.^2));
     cov(5,5,a) = sum((rho2) * real(cfac_5.^2));
     cov(1,2,a) = sum((rho2) * real((cfac_1.*cfac_1).*exp(1i.*fourModes.* (stencil(1) - stencil(2)))));
     cov(1,3,a) = sum((rho2) * real((cfac_1.*cfac_3).*exp(1i.*fourModes.* (stencil(1) - stencil(3)))));
     cov(1,4,a) = sum((rho2) * real((cfac_1.*cfac_4).*exp(1i.*fourModes.* (stencil(1)-stencil(4)))));
     cov(1,5,a) = sum((rho2) * real((cfac_1 .* cfac_5).*exp(1i.*fourModes.* (stencil(1)-stencil(5)))));
     cov(2,1,a) = sum((rho2) * real((cfac_2 .* cfac_1).*exp(1i.*fourModes.* (stencil(2) - stencil(1)))));
     cov(2,3,a) = sum((rho2) * real((cfac_2 .* cfac_3).*exp(1i.*fourModes.* (stencil(2) - stencil(3)))));
     cov(2,4,a) = sum((rho2) * real((cfac_2 .* cfac_4).*exp(1i.*fourModes.* (stencil(2)-stencil(4)))));
     cov(2,5,a) = sum((rho2) * real((cfac_2 .* cfac_5).*exp(1i.*fourModes.* (stencil(2)-stencil(5)))));
     cov(3,1,a) = sum((rho2) * real((cfac_3 .* cfac_1).*exp(1i.*fourModes.* (stencil(3) - stencil(1)))));
     cov(3,2,a) = sum((rho2) * real((cfac_3 .* cfac_2).*exp(1i.*fourModes.* (stencil(3) - stencil(2)))));
     cov(3,4,a) = sum((rho2) * real((cfac_3 .* cfac_4).*exp(1i.*fourModes.* (stencil(3)-stencil(4)))));
     cov(3,5,a) = sum((rho2) * real((cfac_3 .* cfac_5).*exp(1i.*fourModes.* (stencil(3)-stencil(5)))));
     cov(4,1,a) = sum((rho2) * real((cfac_4 .* cfac_1).*exp(1i.*fourModes.* (stencil(4)-stencil(1)))));
     cov(4,2,a) = sum((rho2) * real((cfac_4 .* cfac_2).*exp(1i.*fourModes.* (stencil(4)-stencil(2)))));
     cov(4,3,a) = sum((rho2) * real((cfac_4 .* cfac_3).*exp(1i.*fourModes.* (stencil(4)-stencil(3)))));
     cov(4,5,a) = sum((rho2) * real((cfac_4 .* cfac_5).*exp(1i.*fourModes.* (stencil(4)-stencil(5)))));
     cov(5,1,a) = sum((rho2) * real((cfac_5 .* cfac_1).*exp(1i.*fourModes.* (stencil(5)-stencil(1)))));
     cov(5,2,a) = sum((rho2) * real((cfac_5 .* cfac_2).*exp(1i.*fourModes.* (stencil(5)-stencil(2)))));
     cov(5,3,a) = sum((rho2) * real((cfac_5 .* cfac_3).*exp(1i.*fourModes.* (stencil(5)-stencil(3)))));
     cov(5,4,a) = sum((rho2) * real((cfac_5 .* cfac_4).*exp(1i.*fourModes.* (stencil(5)-stencil(4)))));


     %Define D 
    D = M'*inv(cov(:,:,a))*M;
    
    %Calculate gamma_dot. This is the gamma_dot that will ensure operation
    %at a desired false alarm rate. From Aditya's thesis. We need to
    %normalize the result from Qinv by using the standard deviation and
    %mean. 
%      nonstandard_result = real(Qinv(pfa));
%      standard_result = nonstandard_result*D ;
%      gamma_dot(a) = standard_result * sqrt(D);
%     
    gamma_dot(a) = Qinv(pfa) * sqrt(D);
    %define Y
    Y = zeros(5,1);
    jmpFncCfs_1 = 1i * fHat_alt' .* sign(fourModes)' .* cfac_1';
    jmpFncCfs_2 = 1i * fHat_alt' .* sign(fourModes)' .* cfac_2';
    jmpFncCfs_3 = 1i * fHat_alt' .* sign(fourModes)' .* cfac_3';
    jmpFncCfs_4 = 1i * fHat_alt' .* sign(fourModes)' .* cfac_4';
    jmpFncCfs_5 = 1i * fHat_alt' .* sign(fourModes)' .* cfac_5';
    Y_1 = real(Exponent_factor_1 * jmpFncCfs_1');
    Y_2 = real(Exponent_factor_2 * jmpFncCfs_2');
    Y_3 = real(Exponent_factor_3 * jmpFncCfs_3');
    Y_4 = real(Exponent_factor_4 * jmpFncCfs_4');
    Y_5 = real(Exponent_factor_5 * jmpFncCfs_5');
    Y(1) = Y_1; Y(2) = Y_2; Y(3) = Y_3; Y(4) = Y_4; Y(5) = Y_5; 
    Y = Y;
    %Calculate test statistic
    test_stat(a) = (M'*inv(cov(:,:,a))*Y)';
    %compare test statistic to gamma to detect jump
    % Define the noisy template approximation.
    
end

    Noise_jmpFncCfs = 1i*fHat' .* sign(fourModes)' .* cfac_t';
    Noise_template =  fourMat*Noise_jmpFncCfs' ;
    jump_candidate = NaN(length(FineX),1);
    for a = 1:length(FineX)
        if abs(test_stat(a)) > abs(real(gamma_dot(a)))
        jump_candidate(a) = 1;
        end
    end
%     figure;
%     plot(FineX,jump_candidate); grid on;
    
for a = 7:length(FineX)-6
    if jump_candidate(a) == 1
      if test_stat(a) > test_stat(a+5)
       if (test_stat(a)) > (test_stat(a+4))
            if (test_stat(a)) > (test_stat(a+3))
            if (test_stat(a)) > (test_stat(a+2))
                if (test_stat(a)) > (test_stat(a+1))
                    if (test_stat(a)) > (test_stat(a-1))
                        if (test_stat(a)) > (test_stat(a-2))
                            if (test_stat(a)) > (test_stat(a-3)); 
                                if (test_stat(a)) > (test_stat(a-4))
                                    if test_stat(a) > test_stat(a-5)
                                        jump_detect(a) = -real(template(a));
                                        test_stat_adv(a) = abs(test_stat(a)) - gamma_dot(a);
                                    end
                                end
                            end
                            end
                            end
                        end
                    end
            end
       end
      end
    end
end
     

    
    if jump_candidate(a) ~=1
    jump_detect(a) = NaN;
    end

for a = 4:length(FineX)-3;
    if  abs(jump_detect(a)) < .05
        jump_detect(a) = NaN;
    end
end

jump_detect(1) = NaN;
jump_detect(2) = NaN;
jump_detect(3) = NaN;
jump_detect(4) = NaN;
jump_detect(length(FineX)) = NaN;
jump_detect(length(FineX)-1) = NaN;
jump_detect(length(FineX)-2) = NaN;
jump_detect(length(FineX)-3) = NaN;


figure;
plot(FineX, fx, 'k', FineX, jump_detect, 'ro'); title('Detection of jumps for sawtooth function'); legend('f(x)', 'jump');
figure; 
plot(FineX,test_stat_adv); grid on; 
title('Difference in Test Stat')


%% Now run the second threshold test. We take the x-values 
%  found from the test above and compare their template values 
%  against a scaled threshold that will give us a final set 
%  of jump locations and values.

% ThirdThresholdValues = NaN(length(x),1);
% Second_Threshold = prctile(test_stat_adv,98.5); 
% Final_Jump_Values = NaN(length(FineX),1);
% for a = 1:length(FineX)
%    if jump_detect(a) ~= NaN
%        if abs(test_stat_adv(a)) > Second_Threshold
%            Final_Jump_Values(a) = -template(a);
%            ThirdThresholdValues(a) = 1; 
%            
%        end   
%    else 
%         Final_Jump_Values = 0;
%    end
% end
% figure;
% plot(FineX,fx,FineX,Final_Jump_Values,'ko'); title('Predicted Jumps')
% 
% 
% % Apply the second threshold test
% Second_Threshold = .8; 
% Final_Jump_Values = NaN(length(FineX),1);
% ThirdThresholdValues = NaN(length(FineX),1);
% 
% for a = 1:length(FineX)
%    if jump_detect(a) ~= NaN
%        if abs(jump_detect(a)) > Second_Threshold
%            Final_Jump_Values(a) = template(a);
%            ThirdThresholdValues(a) = 1;
%            
%        end   
%    else 
%         Final_Jump_Values = 0;
%    end
% end
% figure; plot(FineX,fx,FineX,Final_Jump_Values, 'ko'); grid on;
% title('Final Predicted Values') 
% legend('f(x)','Predicted Jumps')
% 
% %% %% We now want to apply the third threshold test, which will use a minmod function to 
% % to rule out the points whose signs based on the jump responses don't
% % match.
% 
% % First calculate the jump response with the exponential. 
% % Template or jump response
% jmpFncCfs_exp =  1i*fHat' .* sign(fourModes)' .* cfe';
% template_exp =  fourMat*jmpFncCfs_exp' ;
% figure; 
% plot(FineX,real( template_exp)); grid on; title('Jump Approximation, Exponential')
% 
% % Then with the trigonometric concentration factor
% 
% jmpFncCfs_trig =  1i*fHat' .* sign(fourModes)' .* cfac_t';
% template_trig =  fourMat*jmpFncCfs_trig' ;
% figure; 
% plot(FineX,real( template_trig)); grid on; title('Jump Approximation, Polynomial')
% 
% 
% % Then with the new concentration factor.
% jmpFncCfs_CFC =  1i*fHat' .* sign(fourModes)' .* sig';
% template_CFC =  fourMat*jmpFncCfs_CFC' ;
% figure; 
% plot(FineX,real( template_CFC)); grid on; title('Jump Approximation, New Factor')
% hold on
% 
% % Then with the first-order polynomial concentration factor.
% jmpFncCfs_poly =  1i*fHat' .* sign(fourModes)' .* cfp';
% template_poly =  fourMat*jmpFncCfs_poly' ;
% figure; 
% plot(FineX,real( template_poly)); grid on; title('Jump Approximation, Polynomial')
% hold off;
% 
% PostFinalValues = NaN(length(FineX),1);
% for a =1:length(FineX)
%     if ThirdThresholdValues(a) == 1;
%         if sign(real(template_trig(a))) >0 && sign(real(template_poly(a))) >0 && sign(template_CFC(a)) >0;
%                 %% Put in the trig once the signs are the same.
%                 PostFinalValues(a) = -template_trig(a);
% %                 Comparison(1) = template_trig(a); 
% %                 Comparison(2) = template_poly(a); 
% %                 Comparison(3) = template_CFC(a);
% %                 PostFinalValues(a) = -min(Comparison);   
%         end
%         if sign(real(template_trig(a))) <0 && sign(real(template_poly(a))) <0 && sign(template_CFC(a)) <0;
%             PostFinalValues(a) = -template_trig(a);
% %                 Comparison(1) = template_trig(a); 
% %                 Comparison(2) = template_poly(a); 
% %                 Comparison(3) = template_CFC(a);
% %                 PostFinalValues(a) = -max(Comparison);   
%         end
%     end
% end
% 
% 
% figure;
% plot(FineX,fx,FineX,PostFinalValues,'ko'); title('Final Predicted Values')
% legend('fx','Predicted Jumps')