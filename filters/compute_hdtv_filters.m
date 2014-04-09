function G = compute_hdtv_filters(degree)
% Compute filters cooresponding to discrete derivative operators
% defined by the tensor product of derivative of B-spline operators
% Results of output for degree = 1:5 are saved as hdtv(n).mat for quick
% loading
    
    if(degree == 1)
        x = [-0.5 0.5];
        d0 = BsplineDiff(x,1,0);
        d1 = BsplineDiff(x,1,1);
        
        G = zeros(length(x),length(x),length(x),3);        
        G(:,:,:,1) = tensor_prod(d1,d0,d0);
        G(:,:,:,2) = tensor_prod(d0,d1,d0);
        G(:,:,:,3) = tensor_prod(d0,d0,d1);            
    end    
    
    if(degree == 2)
        x = [-1 0 1];
        d0 = BsplineDiff(x,2,0);
        d1 = BsplineDiff(x,2,1);
        d2 = BsplineDiff(x,2,2);
        
        G = zeros(length(x),length(x),length(x),6);        
        G(:,:,:,1) = tensor_prod(d2,d0,d0);
        G(:,:,:,2) = tensor_prod(d0,d2,d0);
        G(:,:,:,3) = tensor_prod(d0,d0,d2);
        G(:,:,:,4) = tensor_prod(d1,d1,d0);
        G(:,:,:,5) = tensor_prod(d1,d0,d1);
        G(:,:,:,6) = tensor_prod(d0,d1,d1);               
    end
    
    if(degree == 3)
        x = [-1.5 -0.5 0.5 1.5];        
        d0 = BsplineDiff(x,3,0);
        d1 = BsplineDiff(x,3,1);
        d2 = BsplineDiff(x,3,2);
        d3 = BsplineDiff(x,3,3);        
        
        G = zeros(length(x),length(x),length(x),10);        
        G(:,:,:,1) = tensor_prod(d3,d0,d0);
        G(:,:,:,2) = tensor_prod(d0,d3,d0);
        G(:,:,:,3) = tensor_prod(d0,d0,d3);
        G(:,:,:,4) = tensor_prod(d2,d1,d0);
        G(:,:,:,5) = tensor_prod(d2,d0,d1);
        G(:,:,:,6) = tensor_prod(d1,d2,d0);
        G(:,:,:,7) = tensor_prod(d0,d2,d1);
        G(:,:,:,8) = tensor_prod(d1,d0,d2);
        G(:,:,:,9) = tensor_prod(d0,d1,d2);
        G(:,:,:,10) = tensor_prod(d1,d1,d1);
    end    
    
end