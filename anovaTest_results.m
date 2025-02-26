function anovaTest_results(M)
    p = 0;
    disp('One way ANOVA test results ...');
    for c = 1:size(M, 1)
        if M(c, end) <= 0.005
            disp([num2str(M(c, 1)) ' - ' num2str(M(c, 2)) ': Corrected p-value < 0.005']);
            p = p + 1;
        elseif M(c, end) <= 0.05
            disp([num2str(M(c, 1)) ' - ' num2str(M(c, 2)) ': Corrected p-value < 0.05']);
            p = p + 1;
        end
    end
    if p == 0
        disp('Non of the combinations has significant difference ...');
    end