function interactiveConductance3Dslice(fa, dim_size, potentials, mask_index, dim, slice)
% Interactive conductance map from 3D data, computing one slice
%
% Author: Aina Frau-Pascual


% Conductance map
if dim==3
    fa2 = squeeze(fa(:,:,slice));
    dimx_size = dim_size([1,2]);
elseif dim==2
    fa2 = squeeze(fa(:,slice,:));
    dimx_size = dim_size([1,3]);
else
    fa2 = squeeze(fa(slice,:,:));
    dimx_size = dim_size([2,3]);
end
mnx = prod(dimx_size);

% Interactive conductance map
figure(107); image(fa2','CDataMapping','scaled'); colormap jet; colorbar;
title('FA of the image')

% Interactive conductance map
for i=1:100
    figure(107);
    [x, y] = ginput(1);
    if dim==3
        a1 = round(x); b1 = round(y); c1 = slice;
    elseif dim==2
        a1 = round(x); b1 = slice; c1 = round(y);
    else
        a1 = slice; b1 = round(x); c1 = round(y);
    end
    aux1 = sub2ind(dimx_size, a1, b1);
    p1 = find(mask_index==aux1);
    
    c = zeros(dimx_size);
    for p3=1:mnx
        [a3, b3] = ind2sub(dimx_size, p3);
        if dim==3
            a2 = a3; b2 = b3; c2 = 1;
        elseif dim==2
            a2 = a3; b2 = 1; c2 = b3;
        else
            a2 = 1; b2 = a3; c2 = b3;
        end
        p2 = find(mask_index==p3);
                        
        if not(isempty(p1) || isempty(p2))
            c(a3, b3) = abs(1 / ((potentials(a1,b1,c1,p2) - potentials(a1,b1,c1,p1)) ... 
                             - (potentials(a2,b2,c2,p2) - potentials(a2,b2,c2,p1))));
        end
    end

    % axial
    max_c = max(c(~isinf(c)));
    figure(105); imagesc((c)'); caxis([0 max_c]); colormap(jet); colorbar;
    title('Conductance map for a certain sink point')
end
