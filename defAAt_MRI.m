%Defines Fourier Undersampling operators A, At for Compressed Sensing MRI
%Inputs: acc = acceleration factor, siz = image size
%Outputs: A = measurement operator, At = measurement operator transpose
function [A At S]= defAAt_MRI(options)

    acc = options.acc;
    siz = options.siz;
    density = options.density;
    
    mask3d = zeros(siz);
    if(density == 0)     %uniform random sampling
        mask3d = randSampling(options.siz(1),options.siz(2),options.siz(3),1/acc); 
        mask3d(1,1,1)=1;  %include dc component
    elseif(density == 1)  %density compensation in x-y plane
        for i=1:siz(3)
            pdf = genPDF(siz,5,1/acc);
            mask = im2double(genSampling(pdf,2,5)); 
            mask3d(:,:,i) = mask;
        end
    end
    
    S = find(mask3d~=0);
   
    A = @(z) A_fhp(z,S,siz);
    At = @(z) At_fhp(z,S,siz);

    function A = A_fhp(z, S, siz)
        A = zeros(siz,'double');
        p = 1/sqrt(siz(1)*siz(2)*siz(3))*fftn(z);
        A(S) = p(S);
        A = A(:);
    end

    function At = At_fhp(z, S, siz)
        z = reshape(z,siz);
        p = zeros(siz,'double');
        p(S) = z(S);
        At = sqrt(siz(1)*siz(2)*siz(3))*ifftn(p);
    end

end