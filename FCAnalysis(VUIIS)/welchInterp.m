%% function to interpolate welch power plot
function sZ = welchInterp(l, fs)
    sZ = 0.000001;
    while (length(0:sZ:fs) >= l)
        sZ = sZ + 0.000001;
    end
    sZ = sZ - 0.000001;