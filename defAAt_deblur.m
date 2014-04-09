%Defines Fourier Undersampling operators A, At for Compressed Sensing MRI
%Inputs: acc = acceleration factor, siz = image size
%Outputs: A = measurement operator, At = measurement operator transpose
function [A At]= defAAt_deblur(siz)

    %Define gaussian filter
    hsiz = [5,5,5];
    sig = [1,1,1]; %hsiz/(4*sqrt(2*log(2)));
    hsiz = (hsiz-1)/2;
    [x,y,z] = ndgrid(-hsiz(1):hsiz(1),-hsiz(2):hsiz(2),-hsiz(3):hsiz(3));
    h = exp(-(x.*x/2/sig(1)^2 + y.*y/2/sig(2)^2 + z.*z/2/sig(3)^2));
    h = h/sum(h(:));
    H = psf2otf(h,siz);
    
    A = @(x) ifftn(fftn(x).*H);
    At = @(x) ifftn(fftn(x).*conj(H));
end