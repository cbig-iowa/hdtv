function [ D, Dt ] = defDDt_generic(G,siz)
%DEFDDT_GENERIC Return function handles for generic derivative operators
%loaded as G = 4-D array containing many 3-D filters, siz is size of 3-D
%image to be filtered.

    D = @(X) DDef(X); 
    Dt = @(Y) DtDef(Y); 
    
    %Center and pad filters
    [wx wy wz n] = size(G);
    f = zeros([siz, n]);
    indx = (1:wx)-ceil((wx+1)/2); indx = mod(indx,siz(1))+1; 
    indy = (1:wy)-ceil((wy+1)/2); indy = mod(indy,siz(2))+1;
    indz = (1:wz)-ceil((wz+1)/2); indz = mod(indz,siz(3))+1;    
    f(indx,indy,indz,:) =  G; 
    F = fftn(f);

    %Compute D
    function DX = DDef(X)
        DX = zeros([siz, n],'double');
        for i = 1:n
            DX(:,:,:,i) = ifftn(F(:,:,:,i).*fftn(X));
        end
    end

    %Compute adjoint of D
    function Z = DtDef(Y)
        Z = zeros(siz,'double');
        for i = 1:n
            Z = Z + ifftn(conj(F(:,:,:,i)).*fftn(Y(:,:,:,i)));
        end
    end

end