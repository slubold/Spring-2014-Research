
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Enter number of Fourier modes and grid points.

N = 64;
k = (-N:N)';
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
C = (pi)./(integral(normal,(1/N),1-(1/N)));
Concentration_exp = @(x) (C*x*exp(1./(alpha.*x.*(x-1))));

cfe = zeros(2*N+1,1);    
cfp = zeros(2*N+1,1);
cft = zeros(2*N+1,1);
for i = 1:2*N+1
   s = abs(k(i))./(N);
   cfe(i) = Concentration_exp(s);
   cfp(i) = Concentration_poly(s);
end

%We use the trigonometric/Gibbs concentration factor
Sipi = 1.85193705198247;					% Value of Si(pi), the sine 
                                            % integral                                            
cfac_t = pi * sin(pi*abs(k)/N)/Sipi; 

% % Set cfe(n/2) and cfe(-n/2) = 0
cfe(1) = 0;
cfe(length(cfe)) = 0;

figure(1);
plot(k,cfe,'r',k,cfp,'g',k,cfac_t,'b');
title('Concentration Factors'); grid on;
legend('Exponential', 'Polynomial', 'Trigonometric');

sig = UpdatedConcDesign;
cfac_1 = cfac_t;
cfac_2 = cfac_t;
cfac_3 = sig;
cfac_4 = cfac_t;
cfac_5 = sig;

%% Now we enter the function we want to run the jump detection on, and 
%  calculate the Fourier coefficients.

fx = ((-x-pi).*(x < 0) + (-x+pi).*(x>=0))/(2*pi); fx(x==0) = NaN;
figure; plot(x,fx,'k'); grid on;

for a = 1:2*N+1
    fhat(a) = 1/(2*pi*k(a)*1i); 
end

fhat(N+1) = 0;
figure; plot(k,fhat,'k'); grid on; title('Fourier Coefficients'); 
Func_Approx = fhat * (exp(-1i*x*k'))';
figure; plot(x,imag(Func_Approx),'k'); grid on;
title('Fourier Approximation'); xlabel('x'); ylabel('y');

% Now create the jump repsonse:
fourMat = exp(1i*x*k.' );
jmpFncCfs =  1i*fhat .* sign(k)' .* cfe';
template =  fourMat*jmpFncCfs' ;
figure; 
plot(x,real(template)); grid on; title('Jump Response');
xlabel('x'); ylabel('y');

%% We now turn to the detection of jumps.

%% Noise characteristics
% We consider zero mean, additive white Gaussian nose of variance rho^2
rho2 = 1/N^2;
n = sqrt(rho2)*randn(2*N+1,1) + 1i*sqrt(rho2)*randn(2*N+1,1);
fHat_alt  = fhat' + n;

%delta
delta = 1/(2*N + 1);
%define a probability for false alarms
pfa = .01;

jump_detect = zeros(length(x),1);
test_stat_adv = zeros(length(x),1);
gamma_dot = NaN(length(x),1);
test_stat = NaN(length(x),1);
cov = zeros(5,5,length(x));
for a = 1:length(x)

% First calculate the nonshifted jump responses:  
    %Define stencil
    stencil(1) = x(a) - 2*delta;
    stencil(2) = x(a) - delta;
    stencil(3) = x(a);
    stencil(4) = x(a) + delta;
    stencil(5) = x(a) + 2*delta;
    %Define M
     M = zeros(5,1);
     Exponent_factor_1 = exp( 1i*(stencil(1)) * k.');
     Exponent_factor_2 = exp( 1i*(stencil(2)) * k.');
     Exponent_factor_3 = exp( 1i*(stencil(3)) * k.');
     Exponent_factor_4 = exp( 1i*(stencil(4)) * k.');
     Exponent_factor_5 = exp( 1i*(stencil(5)) * k.');
     jmpFncCfs_1 = 1i * fhat .* sign(k)' .* cfac_1';
     jmpFncCfs_2 = 1i * fhat .* sign(k)' .* cfac_2';
     jmpFncCfs_3 = 1i * fhat .* sign(k)' .* cfac_3';
     jmpFncCfs_4 = 1i * fhat .* sign(k)' .* cfac_4';
     jmpFncCfs_5 = 1i * fhat .* sign(k)' .* cfac_5';
        M_1 = real(Exponent_factor_1 * jmpFncCfs_1');
        M_2 = real(Exponent_factor_2 * jmpFncCfs_2');
        M_3 = real(Exponent_factor_3 * jmpFncCfs_3');
        M_4 = real(Exponent_factor_4 * jmpFncCfs_4');
        M_5 = real(Exponent_factor_5 * jmpFncCfs_5');
        M(1) = M_1; M(2) = M_2; M(3) = M_3; 
        M(4) = M_4; M(5) = M_5;
    %Define covariance matrix

     cov(1,1,a) = sum((rho2) * real(cfac_1.^2));
     cov(2,2,a) = sum((rho2) * real(cfac_2.^2));
     cov(3,3,a) = sum((rho2) * real(cfac_3.^2));
     cov(4,4,a) = sum((rho2) * real(cfac_4.^2));
     cov(5,5,a) = sum((rho2) * real(cfac_5.^2));
     cov(1,2,a) = sum((rho2) * real((cfac_1 .*cfac_1).*exp(1i.*k.* (stencil(1) - stencil(2)))));
     cov(1,3,a) = sum((rho2) * real((cfac_1 .*cfac_3).*exp(1i.*k.* (stencil(1) - stencil(3)))));
     cov(1,4,a) = sum((rho2) * real((cfac_1 .*cfac_4).*exp(1i.*k.* (stencil(1)-  stencil(4)))));
     cov(1,5,a) = sum((rho2) * real((cfac_1 .* cfac_5).*exp(1i.*k.* (stencil(1)-stencil(5)))));
     cov(2,1,a) = sum((rho2) * real((cfac_2 .* cfac_1).*exp(1i.*k.* (stencil(2) - stencil(1)))));
     cov(2,3,a) = sum((rho2) * real((cfac_2 .* cfac_3).*exp(1i.*k.* (stencil(2) - stencil(3)))));
     cov(2,4,a) = sum((rho2) * real((cfac_2 .* cfac_4).*exp(1i.*k.* (stencil(2)- stencil(4)))));
     cov(2,5,a) = sum((rho2) * real((cfac_2 .* cfac_5).*exp(1i.*k.* (stencil(2)-stencil(5)))));
     cov(3,1,a) = sum((rho2) * real((cfac_3 .* cfac_1).*exp(1i.*k.* (stencil(3) - stencil(1)))));
     cov(3,2,a) = sum((rho2) * real((cfac_3 .* cfac_2).*exp(1i.*k.* (stencil(3) - stencil(2)))));
     cov(3,4,a) = sum((rho2) * real((cfac_3 .* cfac_4).*exp(1i.*k.* (stencil(3)-stencil(4)))));
     cov(3,5,a) = sum((rho2) * real((cfac_3 .* cfac_5).*exp(1i.*k.* (stencil(3)-stencil(5)))));
     cov(4,1,a) = sum((rho2) * real((cfac_4 .* cfac_1).*exp(1i.*k.* (stencil(4)-stencil(1)))));
     cov(4,2,a) = sum((rho2) * real((cfac_4 .* cfac_2).*exp(1i.*k.* (stencil(4)-stencil(2)))));
     cov(4,3,a) = sum((rho2) * real((cfac_4 .* cfac_3).*exp(1i.*k.* (stencil(4)-stencil(3)))));
     cov(4,5,a) = sum((rho2) * real((cfac_4 .* cfac_5).*exp(1i.*k.* (stencil(4)-stencil(5)))));
     cov(5,1,a) = sum((rho2) * real((cfac_5 .* cfac_1).*exp(1i.*k.* (stencil(5)-stencil(1)))));
     cov(5,2,a) = sum((rho2) * real((cfac_5 .* cfac_2).*exp(1i.*k.* (stencil(5)-stencil(2)))));
     cov(5,3,a) = sum((rho2) * real((cfac_5 .* cfac_3).*exp(1i.*k.* (stencil(5)-stencil(3)))));
     cov(5,4,a) = sum((rho2) * real((cfac_5 .* cfac_4).*exp(1i.*k.* (stencil(5)-stencil(4)))));


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
    jmpFncCfs_1 = 1i * fHat_alt' .* sign(k)' .* cfac_1';
    jmpFncCfs_2 = 1i * fHat_alt' .* sign(k)' .* cfac_2';
    jmpFncCfs_3 = 1i * fHat_alt' .* sign(k)' .* cfac_3';
    jmpFncCfs_4 = 1i * fHat_alt' .* sign(k)' .* cfac_4';
    jmpFncCfs_5 = 1i * fHat_alt' .* sign(k)' .* cfac_5';
    Y_1 = real(Exponent_factor_1 * jmpFncCfs_1');
    Y_2 = real(Exponent_factor_2 * jmpFncCfs_2');
    Y_3 = real(Exponent_factor_3 * jmpFncCfs_3');
    Y_4 = real(Exponent_factor_4 * jmpFncCfs_4');
    Y_5 = real(Exponent_factor_5 * jmpFncCfs_5');
    Y(1) = Y_1; Y(2) = Y_2; Y(3) = Y_3; Y(4) = Y_4; Y(5) = Y_5; 
 
    %Calculate test statistic
    test_stat(a) = (M'*inv(cov(:,:,a))*Y)';
    %compare test statistic to gamma to detect jump
    % Define the noisy template approximation.
    
end

    Noise_jmpFncCfs = 1i*fhat .* sign(k)' .* cfac_t';
    Noise_template =  fourMat*Noise_jmpFncCfs';
    jump_candidate = NaN(length(x),1);
    for a = 1:length(x)
           if abs(test_stat(a)) > abs(real(gamma_dot(a)))
           jump_detect(a) = Noise_template(a);
           end
    end
        
for a = 1:length(x)
    if  jump_detect(a) ==  0
        jump_detect(a) = NaN;
    end
end
figure;
plot(x, fx, 'k', x, jump_detect, 'ro'); 
title('Detection of jumps for sawtooth function'); 
legend('f(x)', 'jump');

for a = 2:length(x)-1
    if jump_detect(a) == NaN
    jump_detect2(a) = NaN;
    else 
            if abs(test_stat(a)) > abs(test_stat(a+1))
                if abs(test_stat(a)) > abs(test_stat(a-1))
                    jump_detect2(a) = Noise_template(a);
                else
                    jump_detect2(a) = NaN;
            
                end
            else 
                   jump_detect2(a) = NaN;
            end
    end
end
jump_detect2(1) =0;
jump_detect2(length(x)) = 0;
figure;
plot(x, fx, 'k', x, jump_detect2, 'ro'); grid on;
title('Detection of jumps for sawtooth function'); 
legend('f(x)', 'jump');   

%     Noise_jmpFncCfs = 1i*fhat .* sign(k)' .* cfe';
%     Noise_template =  fourMat*Noise_jmpFncCfs' ;
%     jump_candidate = NaN(length(x),1);
%     for a = 1:length(x)
%         if abs(test_stat(a)) > abs(real(gamma_dot(a)))
%         jump_candidate(a) = 1;
%         end
%     end
%     figure;
%     plot(FineX,jump_candidate); grid on;
%     
% for a = 7:length(x)-6
%     if jump_candidate(a) == 1
%         if test_stat(a) > test_stat(a+5)
%             if (test_stat(a)) > (test_stat(a+4))
%                 if (test_stat(a)) > (test_stat(a+3))
%                     if (test_stat(a)) > (test_stat(a+2))
%                         if (test_stat(a)) > (test_stat(a+1))
%                             if (test_stat(a)) > (test_stat(a-1))
%                                 if (test_stat(a)) > (test_stat(a-2))
%                                     if (test_stat(a)) > (test_stat(a-3)); 
%                                         if (test_stat(a)) > (test_stat(a-4))
%                                             if test_stat(a) > test_stat(a-5)
%                                                 jump_detect(a) = -real(template(a));
%                                                 test_stat_adv(a) = abs(test_stat(a)) - gamma_dot(a);
%                                     end
%                                 end
%                             end
%                             end
%                             end
%                         end
%                     end
%             end
%        end
%        end
%     end
%     end
%      
% 
%     
%     if jump_candidate(a) ~=1
%     jump_detect(a) = NaN;
%     end
% 
%     for a = 4:length(x)-3;
%     if  abs(jump_detect(a)) < .05
%         jump_detect(a) = NaN;
%     end
% end
% 
% jump_detect(1) = NaN;
% jump_detect(2) = NaN;
% jump_detect(3) = NaN;
% jump_detect(4) = NaN;
% jump_detect(length(x)) = NaN;
% jump_detect(length(x)-1) = NaN;
% jump_detect(length(x)-2) = NaN;
% jump_detect(length(x)-3) = NaN;

