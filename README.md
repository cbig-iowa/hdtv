HDTV3D
====
This is a MATLAB implementation of Higher Degree Total Variation (HDTV) regularization penalty for use in 3-D image denoising, deblurring, compressed sensing MRI recovery, and other inverse problems. HDTV is an extension of the popular total variation (TV) penalty to higher degree derivatives. Whereas TV promotes images with sparse gradients, HDTV regularization promotes the sparsity of all *n*th degee directional derivatives in the image.

For more information, see: (https://research.engineering.uiowa.edu/cbig/content/generalized-hdtv)

##Code Basics
To begin, open the MATLAB script "start.m". Running this script unaltered will demonstrate the use of HDTV in denoising a 256x256x16 axial brain MR image. To adjust the amount of regularization, change the variable `lambda`; increasing `lambda` will increase the amount of regularization (i.e smooth more), decreasing `lambda` will apply less regularization. *In all settings, the `lambda` parameter must be tuned to obtain optimal results*.

##Changing the Mode
This package is set-up to perform 3-D image denoising, deblurring, and compressed sensing MRI recovery.

###Denoising
The denoising mode is enabled by default. To change the level of noise, change the `SNR` variable. This will be the approximate signal-to-noise ratio of the noisy image.

###Deblurring
To enable the deblurring mode, uncomment the following lines:
    %deblurring
    options.complex = 0;
    [A At] = defAAt_deblur(options.siz);
By default the blurring kernel is an isotropic gaussian filter of size 5x5x5 with standard deviation 1. To change the kernel, you may alter the function defintion in "defAAt_deblur.m". Note: by default noise is also added to the blurred image.

###CS-MRI Recovery
To enable CS-MRI mode, uncomment the following lines:
    %compressed sensing MRI/fourier undersamping setting
    options.acc = 1.65; %acceleration factor
    options.complex = 1; %flag measurement operator as complex-valued
    options.density = 1; %use density compensated random samples (=1) 
                        %or uniform random samples (=0)
    [A At S] = defAAt_MRI(options); %S = k-space undersampling mask 
`options.acc` sets the acceleration factor (=1/(undersampling amount)). Setting `options.complex = 1` ensures that the noise added to the k-space measurements is complex. The flag `options.density` determines the k-space undersampling pattern; see the "defAAt_MRI.m" for more details. The array S gives the k-space sampling locations, where (1,1,1) is the center of k-space.

##Optimization Parameters
Aside from the regulization parameter `lambda`, the optimization parameters in the struct `options` should not need tuning.
To change the degree of derivative used in the HDTV penalty, change `options.degree`. Currently, only degrees 1,2, and 3 are supported. The default (2) should be appropriate in most cases. 

If speed is a concern, the execution time of the algorithm can be reduced by decreasing `options.nsamples`--the number of samples used in the spherical quadrature scheme--or by decreasing `options.Nouter` and `options.Ninner`--the number of outer and inner iterations used in our alternating minimization scheme. However, this could result in small losses in SNR.

##References
1. Y. Hu, G. Ongie, S. Ramani, M. Jacob, "Generalized Higher Degree Total Variation", IEEE Trans. Image Proccessing, 2014 (in press).
2. G. Ongie, Y. Hu, M. Jacob, "Higher degree total variation for 3-D  image recovery", ISBI, Beijing, China, 2014.
3. Y. Hu, M. Jacob,"Higher degree total variation (HDTV) regularization for image recovery", IEEE Trans. Image Processing, pp 2559-2571, Vol 21, No 5, May 2012.