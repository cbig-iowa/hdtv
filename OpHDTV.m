function [X, earray] = OpHDTV(b,A,At,lambda,options)
% OPHDTV: Solves HDTV regularized inverse problems with an alternating 
% minimization algorithm. Returns X, the recovered image, and earray,
% the values of the cost function for all iterations.
%
% Set degree of HDTV operator with options.degree
% ex: degree = 1 ==> equivalent to TV
% ex: degree = 2 ==> HDTV2 penalty using 2nd degree directional derivatives

Nouter = options.Nouter;
Ninner = options.Ninner;
beta = options.beta;
bfactor = options.bfactor;
siz = options.siz;
tol = options.tol;
degree = options.degree;
nsamples = options.nsamples;

%Define AtA
p_image = zeros(siz,'double'); p_image(1,1,1) = 1;
AtA = fftn(At(A(p_image)));

%Define derivative operators and steering functions
load(['filters/hdtv',num2str(degree),'.mat']); %loads filter array as G
[D Dt] = defDDt_generic(G,siz); %build derivatve ops from filters G
nfilters = size(G,4);

[pts weights] = getLebedevQuad(nsamples); %get quadrature samples/weights
su = steer(pts,degree); 
su = su.*repmat(weights,[nfilters,1]); %multiply by quadrature weights
                                       %note: weights must be positive
%Precompute quantites
C = su*su.';

Dvec = reshape(D(p_image), [siz(1)*siz(2)*siz(3), nfilters]);
CD = reshape(Dvec*C,[siz,nfilters]);
DtCD = fftn(Dt(CD));
Atb = At(b);
X = Atb; %Initialize X

% Begin alternating minimization alg.
earray = zeros(1,Nouter*Ninner);
ind = 0;
for i=1:Nouter 
    for ii=1:Ninner        
        %Calculate error        
        DX = D(X); %4-d array of all partial derivatives
        DXvec = reshape(DX, [siz(1)*siz(2)*siz(3), nfilters]);
        DD = DXvec*su; %directional derivatives (#pixels x #directions)
        
        diff = A(X)-b;
        error = norm(diff(:)).^2 + lambda*sum(abs(DD(:))); %objective function
        ind = ind+1;
        earray(ind) = error;

        %Convergence test
        if ii > 2
            if abs((earray(ind)-earray(ind-1)))<tol
                break;
            end
        end          
                
        %Shrinkage step
        shrinkDD = shrink(DD,1/beta); %shrinkage of derivative mag.
        DD = shrinkDD.*DD; 
        DX = reshape(DD*su.',[siz, nfilters]); %compute projected shrinkage
        
        %Inversion step 
        F1 = fftn(2*Atb + lambda*beta*Dt(DX));
        F2 = lambda*beta*DtCD + 2*AtA;
        X = ifftn(F1./F2);
        
    end
    
    %Increase continuation parameter
    beta = beta*bfactor;
    
    fprintf('Iteration %d of %d complete.\n',i,Nouter);
end

earray = earray(1:ind); %truncate error array