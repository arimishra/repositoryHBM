function map = contrastAdjust(map)
    map = map/max(map(:));
    thrld = 0.1;
    while (length(find(map > thrld*max(map(:)))) > 7)
        thrld = thrld + 0.1;
    end
    thrld = thrld - 0.1;
    kp = find(map > thrld*max(map(:)));
    map(kp) = thrld*max(map(:)) + 0.16*(1 - thrld)*map(kp);
    end