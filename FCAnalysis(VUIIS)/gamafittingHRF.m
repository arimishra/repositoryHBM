%%{
function [z] = gamafittingHRF(y)
    options = optimset('MaxFunEvals', 40000, 'MaxIter', 40000);
    xlb = [0; 1; 1; 0; 1; 1];
    xub = 40*ones(6, 1);
    x0 = [1; 5; 5; 1; 5; 5];  
    [x, resnorm, res] = lsqcurvefit(@bimodel_twogamma, x0, 1:length(y), y, xlb, xub, options);
    z = bimodel_twogamma(x, 1:length(y));
%}
    %{
    function [z] = gamafittingHRF(y);
        x = 1:length(y);
        p = nlinfit(x, y, @f, [1 40 60]);
        z = f(p, x);
    end
    function y = f(abc, x)
        a = abc(1); b = abc(2); c = abc(3);
        y = c * x.^(a - 1).* exp(-x/b)/(b^a*gamma(a));
    end
    %}