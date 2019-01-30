function interactiveConductance2D(fa, dims, potentials, mask_index)
% Interactive conductance map from 2D data
%
% Author: Aina Frau-Pascual

mnx = prod(dims);

% Conductance map
figure(101); image(fa','CDataMapping','scaled'); colorbar;
title('FA of the image')

% Interactive conductance map
for i=1:100
    figure(101);
    [x, y] = ginput(1);
    a1 = round(x);
    b1 = round(y);
    aux1 = sub2ind(dims, a1, b1);
    p1 = find(mask_index==aux1);
    
    c = zeros(dims);
    for p3=1:mnx
        [a2, b2] = ind2sub(dims, p3);
        p2 = find(mask_index==p3);
        if not(isempty(p1) || isempty(p2))
            c(a2, b2) = abs(1 / ((potentials(a1,b1,p2)-potentials(a1,b1,p1)) ... 
                               - (potentials(a2,b2,p2)-potentials(a2,b2,p1))));
        end
    end

    max_c = max(c(~isinf(c)));
    figure(102); imagesc(c'); caxis([0 max_c]); colormap jet; colorbar;
    title('Conductance map for a certain sink point')
end
