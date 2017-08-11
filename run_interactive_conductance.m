% Conductance map
saveims = false;
figure(101); image(im,'CDataMapping','scaled'); colorbar;
if saveims
    savefig(strcat(name, '.fig'));
    saveas(gcf, strcat(name, '.eps'), 'epsc');
end

% Interactive conductance map
for i=1:100
    figure(101);
    [x, y] = ginput(1);
    a1 = round(y);
    b1 = round(x);
    p1 = sub2ind([Nx, Ny], a1, b1);
    
    c = zeros(Nx, Ny);
    for p2=1:mny
        [a2, b2] = ind2sub([Nx, Ny], p2);
        c(a2, b2) = abs(1 / ((potentials(a1,b1,p2) - potentials(a1,b1,p1)) ... 
                         - (potentials(a2,b2,p2) - potentials(a2,b2,p1))));
    end
    figure(102); imagesc((c), [0,0.5]); colormap jet; colorbar;
    if saveims
        savefig(strcat(name, '_conductance_p', int2str(p1), '.fig'));
        saveas(gcf, strcat(name, '_conductance_p', int2str(p1), '.eps'), 'epsc');
    end
end