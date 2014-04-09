% Solves the HDTV regularized 3-D image recovery problem:
%   min_x ||Ax - b||^2 + lambda*||x||_{HDTV}
% where x is the image to be recovered, A is a measurement operator, 
% b is the vector of (noisy) measurements, lambda is a the regularization 
% parameter and ||-||_{HDTV} is a higher degree total variation semi-norm; 
% see references below for more information.
%
% We give implementations of measurement operators A for compressed 
% sensing MRI, deblurring, and denoising. For comparison purposes, we also 
% include a choice of TV regularization. Note that in all cases tuning 
% of the lambda parameter is necessary to obtain optimal results.
%
% References:
% [1] Y. Hu, G. Ongie, S. Ramani, M. Jacob, "Generalized Higher Degree 
%     Total Variation", IEEE Trans. Image Proccessing, 2014 (in press).
% [2] G. Ongie, Y. Hu, M. Jacob, "Higher degree total variation for 3-D 
%    image recovery", ISBI, Beijing, China, 2014.
% [3] Y. Hu, M. Jacob,"Higher degree total variation (HDTV) regularization
%    for image recovery", IEEE Trans. Image Processing, pp 2559-2571, 
%    Vol 21, No 5, May 2012.
%
% Send Comments/Bug Reports to: gregory-ongie@uiowa.edu
%
% --Greg Ongie, 4/08/2014
%% Load image
load('data/brain.mat');
%x = phantom3d(64);         %Shepp-Logan phantom
%x = phantom3d_linear(64);  %Modified piecewise linear S-L phantom

x0 = x/max(abs(x(:))); %normalize
options.siz = size(x0); 

%% Define measurement operator A and its adjoint At
%compressed sensing MRI/fourier undersamping setting
% options.acc = 1.65; %acceleration factor
% options.complex = 1; %flag measurement operator as complex-valued
% options.density = 1; %use density compensated random samples (=1) 
%                        %or uniform random samples (=0)
% [A At S] = defAAt_MRI(options); %S = k-space undersampling mask 

%deblurring
% options.complex = 0;
% [A At] = defAAt_deblur(options.siz); 

%denoising: A,At are identity operators
options.complex = 0;
A = @(z) z; At = @(z) z; 

%% Define measurements with noise
SNRdB = 20; %set noise level (in terms of SNR dB)
b = A(x0); %measurement vector
stdsignal = sqrt(sum(abs(b(:)).^2)/length(b(:)));
stdnoise = 0.1^(SNRdB/20)*stdsignal;
if options.complex
    %complex noise (for CS MRI)
    bnoise = stdnoise*(randn(size(b)) + 1i*randn(size(b))); 
else
    %real noise (for deblurring/denoising)
    bnoise = stdnoise*randn(size(b)); 
end
b = b + bnoise; 

%% Set optimization parameters
options.Nouter = 7;    %Outer iterations
options.Ninner = 10;   %Inner iterations
options.beta = 20;     %Continuation parameter beta initialization
options.bfactor = 2;   %Multiplication factor for increasing beta
options.tol = 5e-5;    %Inner loop convergence tolerance

%% Run optimization
lambda = 0.003; %Regularization parameter; tune for optimal results

% with TV Regularization (isotropic)
%[x, earray] = OpTV(b,A,At,lambda,options); 

% with HDTV Regularization
options.degree = 2; %degree of HDTV penalty=1,2,or 3
                    %use 2 or 3 for best results
options.nsamples = 86; %number of samples in Lebedev quadrature scheme
                      %choose from: [ 6, 14, 26, 38, 50, 86, 110, 146, 170]
                      %values < 86 may result in suboptimal SNR
tic                  
[x, earray] = OpHDTV(b,A,At,lambda,options);
toc

% calculate SNR
SNR = 10*log10((sum(abs(x0(:)).^2))/(sum(abs(x(:)-x0(:)).^2)));
disp(['Recon SNR = ',num2str(SNR),' dB'])

%% Plots
zslice = ceil(options.siz(3)/2);
figure;
subplot(1,2,1); imshow(abs(x0(:,:,zslice)),[]); title('Original (middle slice)')
subplot(1,2,2); imshow(abs(x(:,:,zslice)),[]); title(['Recon, SNR = ', num2str(SNR)]);

figure;
plot(earray); title('Cost Function vs. Iterations');

