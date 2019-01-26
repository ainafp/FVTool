function conn = computeConductanceMatrix(potentials, masked_atlas)
%%% Compute global connectivity matrix
%
% Author: Aina Frau-Pascual


roi_list0 = unique(masked_atlas);        % consider 0 is background
roi_list = roi_list0(roi_list0>0); 
nroi = size(roi_list,1);
conn = zeros(nroi, nroi);
for p1=1:nroi
    for p2=1:nroi
        if p2>p1
            pot = (potentials(:,:,:,p2) - potentials(:,:,:,p1));
            pot1 = mean(pot(masked_atlas==roi_list(p1)));
            pot2 = mean(pot(masked_atlas==roi_list(p2)));
            conn(p1, p2) = 1 / (pot1 - pot2);
            conn(p2, p1) = 1 / (pot2 - pot1);
        end
    end
end

end

