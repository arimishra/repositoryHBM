function Y_s = post_processMap(Y_s)
    Y_s(isnan(Y_s)) = 0;
    Y_s(isinf(Y_s)) = 0;