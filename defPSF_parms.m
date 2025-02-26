clear all; 
close all;
warning off;
    %% Code written on 07/22/2018 by Arabinda Mishra, PhD 
    %{
    PSF dimensions in both cases are acquired of the t-map/beta maps (spmT_0001/beta_0001/Bstim/PsigC.nii) 
    for stimulus data or correlation/zscore maps (MnCnvt_ar3b/ZMnCnvt_ar3b.nii) using code defPSF_parm.m. 
    This code outputs the major and minor axis of the gaussian fitted point spread finction.
    %}
    pathN = '/Users/mishraa/Desktop/arm/limin_obliqLmn/codePackage/sampleData_rest/';
    cvtMap = 'MnCnvt_ar3b_d2T.nii';
    sub = 4;
    win_widthX = 10;
    win_widthY = 8;
    thrldH = 0.2;
    %% PSF parameter estimation on layer 1-3, resting state and layer 1-2 stimulus data
    if cvtMap(end - 4) == 'T'
        slice = 9;
    elseif cvtMap(end - 4) == 'B'
        slice = 8;
    else
        slice = 7;
    end
    
    load colmap_trP.mat
    Vout = spm_vol([pathN '/anat_img/anatImage.nii']);
    A = spm_read_vols(Vout);
    A = squeeze(A(:, :, slice));
    A = 0.5*A/max(A(:));
    Vout = spm_vol([pathN '/anat_img/' cvtMap]);
    T = spm_read_vols(Vout);
    T = squeeze(T(:, :, slice));
    %% Identifies peak activation
    [i j] = find(T == max(T(:)))
    [x y X Y] = interpAnatFnc(T, A);
    [d1] = voxel_thresholdS(1, (T > thrldH));
    T = d1.*T;
    T = interp2(x, y, T, X, Y, 'nearest');
    d1 = (T > 0);
    
    keeps = find(T > 0);
    mxTmap = max(T(:));
    for c = 1:length(keeps)
        T(keeps(c)) = (T(keeps(c)) - thrldH)/(mxTmap -thrldH);
    end
    dsp = A.*~d1 + 0.52*d1 + 0.48*T;
    imagesc(permute(dsp, [2, 1]));
    colormap(tray);
    axis off;
    x_l = i - 10;
    x_u = i + 10;
    y_l = j - 10;
    y_u = j + 10;
    
    Vout = spm_vol([pathN '/anat_img/anatImage.nii']); 
    map = spm_read_vols(Vout);
    Vout = spm_vol([pathN '/anat_img/' cvtMap]);
    dispIm = spm_read_vols(Vout);
    mul_fact = size(map, 1)/size(dispIm, 1);
    map = squeeze(map(:, :, slice));
    map = map/max(map(:));
    map((x_l - 1)*mul_fact + 1:x_u*mul_fact, (y_l - 1)*mul_fact + 1:(y_l - 1)*mul_fact + 2) = 1;
    map((x_l - 1)*mul_fact + 1:x_u*mul_fact, y_u*mul_fact:y_u*mul_fact + 1) = 1;
    map((x_l - 1)*mul_fact + 1:(x_l - 1)*mul_fact + 2, (y_l - 1)*mul_fact + 1:y_u*mul_fact) = 1;
    map(x_u*mul_fact:x_u*mul_fact + 1, (y_l - 1)*mul_fact + 1:y_u*mul_fact) = 1;
    figure;
    map = permute(map, [2 1]);
    imagesc(map);
    colormap(gray);
    axis off;
    figure;
    disp_Im = permute(squeeze(dispIm(x_l:x_u, y_l:y_u, slice)), [2 1]);
    imagesc(disp_Im);
    [d1] = voxel_thresholdS(1, (disp_Im > thrldH));
    disp_Im = d1.*disp_Im;
    colormap(gray);
    axis off;
    map = map((x_l - 1)*mul_fact + 1:x_u*mul_fact, (y_l - 1)*mul_fact + 1:y_u*mul_fact);
    bigtmap = imresize(disp_Im, [50 50], 'bicubic');
    new_bigtmap_nor = bigtmap./max(bigtmap(:));
    [X Y] = meshgrid(1:size(new_bigtmap_nor, 2), 1:size(new_bigtmap_nor, 1));
    figure;
    mesh(X, Y, flipud(new_bigtmap_nor));
    colormap(jet);
    
    %colorbar;
    xlabel('mm', 'fontname', 'arial', 'fontsize', 24);
    ylabel('mm', 'fontname', 'arial', 'fontsize', 24);
    zlabel('correlation map', 'fontname', 'arial', 'fontsize', 24);
    set(gca,'FontSize', 24);
    set(gcf,'color', 'w');
    thrsl_bigtmap_nor = new_bigtmap_nor;
    thrsl_bigtmap_nor(abs(thrsl_bigtmap_nor) < thrldH) = 0;

    figure;
    contourf(flip(thrsl_bigtmap_nor, 1), [0.3 0.4 0.6]);
    colormap(jet);
    colorbar
    set(gca, 'FontSize', 24);
    set(gcf, 'color', 'w');
    grid on
    rot = imrotate(thrsl_bigtmap_nor, 0, 'bicubic');

    figure;
    imagesc(rot);
    colormap(jet);
    colorbar
    set(gca, 'FontSize', 24);
    set(gcf, 'color', 'w');
    title(['stim/rest map', 'fontsize', 24])

    %% Find the mass point
    figure;
    imagesc(rot);
    colormap(jet);
    colorbar
    set(gca, 'FontSize', 24);
    set(gcf, 'color', 'w');
    %%{
    %% draw ROIs
    r = roipoly;
    M1 = r.*squeeze(rot);
    [centerMass_x, centerMass_y] = find(M1 == max(M1(:)));
    mass_bigtmap_nor = new_bigtmap_nor;
    mass_bigtmap_nor(mass_bigtmap_nor < 0) = 0;
    rot_mass = imrotate(mass_bigtmap_nor, 0, 'bicubic');
    figure;
    imagesc(rot_mass);
    colormap(jet);
    colorbar;
    hold on;
    plot(centerMass_y, centerMass_x, 'w.', 'MarkerSize', 24);
    hold off
    set(gca, 'FontSize', 24);
    set(gcf, 'color', 'w');

    constant_x = centerMass_x - win_widthY + 1:centerMass_x + win_widthY + 1;
    constant_y = centerMass_y - win_widthX + 1:centerMass_y + win_widthX + 1;
    bw_x = rot_mass(centerMass_x, constant_y);
    bw_y = rot_mass(constant_x, centerMass_y)';
    if size(bw_x, 1) > 1
        bw_x = mean(bw_x);
    end
    if size(bw_y, 1) > 1
        bw_y = mean(bw_y);
    end
    figure;
    x = 1:length(bw_x);
    [A, mu, cigma] = mygaussfit(x, bw_x);
    y = A*exp(-((x - mu)./cigma).^2);
    plot(1:length(bw_x), bw_x, '.'); hold on;
    plot(1:length(bw_x), y, 'r'); hold on;
    FWHM_x = 2*sqrt(2*log(2))*cigma*35./512./sqrt(2)
    figure;
    x = 1:length(bw_y);
    [A, mu, cigma] = mygaussfit(x, bw_y);
    y = A*exp(-((x - mu)./cigma).^2);
    plot(1:length(bw_y), bw_y, '.'); hold on;
    plot(1:length(bw_y), y, 'r'); hold on;
    FWHM_y = 2*sqrt(2*log(2))*cigma*35./512./sqrt(2)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
