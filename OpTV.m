function [X, earray] = OpTV(b,A,At,lambda,options)
% OPTV: Solves TV regularized inverse problems with an alternating 
% minimization algorithm. Returns X, the recovered image, and earray,
% the values of the cost function for all iterations.
%
% Based on the fTVd implementation: 
% http://www.caam.rice.edu/~optimization/L1/ftvd/

Nouter = options.Nouter;
Ninner = options.Ninner;
beta = options.beta;
bfactor = options.bfactor;
siz = options.siz;
tol = options.tol;

%Define AtA
p_image = zeros(siz,'double'); p_image(1,1) = 1;
AtA = fft2(At(A(p_image)));

%Define derivative operators D, Dt, and DtD
[D,Dt] = defDDt_TV(siz);
DtD = fft2(Dt(D(p_image)));

%Initialize X
X = At(b);

% Begin alternating minimization alg.
earray = [];
for i=1:Nouter 
    for ii=1:Ninner        
        %Calculate error        
        DX = D(X);
        DX1 = DX(:,:,1); 
        DX2 = DX(:,:,2);
        ADX = sqrt(abs(DX1).^2 + abs(DX2).^2); %isotropic TV           
        diff = A(X)-b;
        error = norm(diff(:)).^2 + lambda*sum(ADX(:)); %objective function
        earray = [earray,error];

        %Convergence test
        if ii > 2
            if abs((earray(end)-earray(end-1)))<tol
                break;
            end
        end          
                
        %Shrinkage step
        shrinkDX = shrink(ADX,1/beta); %shrinkage of gradient mag.
        DX(:,:,1) = shrinkDX.*DX1; 
        DX(:,:,2) = shrinkDX.*DX2;

        %Inversion step 
        F1 = fft2(2*At(b) + lambda*beta*Dt(DX));
        F2 = lambda*beta*DtD + 2*AtA;
        X = ifft2(F1./F2);
        
    end
    
    %Increase continuation parameter
    beta = beta*bfactor;
        
end
