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
for i=1:100
    figure(107);
    [x, y] = ginput(1);
    a1 = round(x);
    b1 = round(y);
    aux1 = sub2ind(dimx_size, a1, b1);
    p1 = find(mask_index==aux1);
    
    c = zeros(dimx_size);
    for p3=1:mnx
        [a2, b2] = ind2sub(dimx_size, p3);
        p2 = find(mask_index==p3);
                        
        if not(isempty(p1) || isempty(p2))
            c(a2, b2) = abs(1 / ((potentials(a1,b1,p2) - potentials(a1,b1,p1)) ... 
                             - (potentials(a2,b2,p2) - potentials(a2,b2,p1))));
        end
    end

    % axial
    figure(105); imagesc((c)'); caxis([0 0.005]); colormap(jet); colorbar;

end
