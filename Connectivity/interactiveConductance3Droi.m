function interactiveConductance3Droi(fa, mask, atlas, dim_size, potentials, dim, slice, max_limit)
% Interactive conductance map
%
% Author: Aina Frau-Pascual


% Conductance map
if dim==3
    fa2 = squeeze(fa(:,:,slice));
    mask0 = squeeze(mask(:,:,slice));
    dimx_size = dim_size([1,2]);
    atlas_slice = squeeze(atlas(:,:,slice));
elseif dim==2
    fa2 = squeeze(fa(:,slice,:));
    mask0 = squeeze(mask(:,slice,:));
    dimx_size = dim_size([1,3]);
    atlas_slice = squeeze(atlas(:,slice,:));
else
    fa2 = squeeze(fa(slice,:,:));
    mask0 = squeeze(mask(slice,:,:));
    dimx_size = dim_size([2,3]);
    atlas_slice = squeeze(atlas(slice,:,:));
end
mnx = prod(dimx_size);

% ROI treatment
roi_list0 = unique(atlas);                       % consider 0 is background
roi_list = roi_list0(roi_list0>0);               % and 2 is WM

% Interactive conductance map
figure(100); image(fa2','CDataMapping','scaled'); colorbar;
for i=1:100
    figure(100);
    [x, y] = ginput(1);
    if dim==3
        a1 = round(x); b1 = round(y); c1 = slice;
    elseif dim==2
        a1 = round(x); b1 = slice; c1 = round(y);
    else
        a1 = slice; b1 = round(x); c1 = round(y);
    end
    roi0 = atlas(a1, b1, c1);
    roi1 = find(roi_list==roi0);
    roi_list
    size(roi_list)
    
    c = zeros(size(atlas(:,:,slice)));
    for roi=1:size(roi_list,1)
        p2_vector = find(atlas_slice==roi_list(roi));
        
        for p2=p2_vector'
            [a3, b3] = ind2sub(dimx_size, p2);
            if dim==3
                a2 = a3; b2 = b3; c2 = slice;
            elseif dim==2
                a2 = a3; b2 = slice; c2 = b3;
            else
                a2 = slice; b2 = a3; c2 = b3;
            end
            %roi
            %size(potentials)
            if not(isempty(roi1) || isempty(roi))
                c(a3, b3) = abs(1 / ((potentials(a1,b1,c1,roi) - potentials(a1,b1,c1,roi1)) ... 
                                   - (potentials(a2,b2,c2,roi) - potentials(a2,b2,c2,roi1))));
            end
        end
    end
    
    jet0 = jet;
    jet0(1,:) = [ 1 1 1 ]; % for white background
    %max_limit = 0.02;
    figure(105); imagesc((c.*mask0)'); colormap(jet); colorbar; caxis([0 max_limit]); 
end

