clear all;
close all;
warning off;
clc
    load area1PSF_R
    digitT = area1T_R;
    digitB = area1B_R;
    digitL = area1L_R;
    
    %% Resting data (FWHM major minor axis)
    figure;
    X = [digitT(:, 1) digitB(:, 1) digitL(:, 1) digitT(:, 2) digitB(:, 2) digitL(:, 2)]; hold on
    boxplot(X), hold on;
    xlabel('Lr1-3 (Mj1-3, Mr1-3, rest)', 'fontname', 'arial', 'fontsize', 16);
    ylabel('FWHM(mm)', 'fontname', 'arial', 'fontsize', 16);
    plot(1:6, mean(X), 'd'), hold on;
    [p, tbl, stats] = anova1(X);
    anovaTest_results(multcompare(stats));
    %manWhitneyTest(X);

    figure;
    %% Resting data Axis ratio
    X = [digitT(:, 3) digitB(:, 3) digitL(:, 3)]; hold on
    boxplot(X), hold on;
    xlabel('Layer 1-3(rest)', 'fontname', 'arial', 'fontsize', 16);
    ylabel('Axis ratio', 'fontname', 'arial', 'fontsize', 16);
    plot(1:3, mean(X), 'd'), hold on;
    [p, tbl, stats] = anova1(X);
    anovaTest_results(multcompare(stats));
    %manWhitneyTest(X);
    
    figure;
    %% Resting data Area enclosed (layer
    X = [digitT(:, 4) digitB(:, 4) digitL(:, 4)]; hold on
    boxplot(X), hold on;
    xlabel('Layer 1-3(rest)', 'fontname', 'arial', 'fontsize', 16);
    ylabel('Area(mm2)', 'fontname', 'arial', 'fontsize', 16);
    plot(1:3, mean(X), 'd'), hold on;
    [p, tbl, stats] = anova1(X);
    anovaTest_results(multcompare(stats));
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    