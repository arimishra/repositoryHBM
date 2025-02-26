%% Program written by Arabinda Mishra, October 20, 2014
close all;
clear all;
warning off;
    %mp_type = 1;        % input 1 - correlated, 2 - anti corelated map and 3 - for both
    %show_slice = 1; %0 - only slice where activation occurs 1 - show all of them
    %thrldL = - 1.0;
    %contrast = 1;
    maps = {'beta_0001', 'spmT_0001', 'PsigC', 'Stim', 'StimH', 'StimT', 'Bstim', 'BstimH', 'BstimT', 'snrRaw', 'Zstim', 'ZstimT'};
    path = '/Users/mishraa/Desktop/arm/limin_obliqLmn/BOLD/HRS_data/';
    %{
    subjs = {'m5793d2_170418_epi21', 'm5793d2_170418_epi23', 'm5793d2_170418_epi26', 'm5793d2_170418_epi31', ...
            'm5793d2_170418_epi32', 'm5793d2_170418_epi34', 'm5793d3_170418_epi22', 'm5793d3_170418_epi24', ...
            'm5793d3_170418_epi32', 'm5793d3_170418_epi35', 'm5793d23_170418_epi20', 'm5793d23_170418_epi25', ...
            'm5793d23_170418_epi28', 'm5793d23_170418_epi33', 'm5793d23_170418_epi36', 'm6224d2_170411_epi27', ...
            'm6224d2_170411_epi30', 'm6224d2_170411_epi32', 'm6224d3_170411_epi26', 'm6224d3_170411_epi29', ...
            'm6224d23_170411_epi28', 'm6224d23_170411_epi31', 'm6268d2_170419_epi19', 'm6268d2_170419_epi21', ...
            'm6268d3_170419_epi20', 'm6268d3_170419_epi22', 'm6268d23_170419_epi17', 'm6268d23_170419_epi18', ...
            'm6268d23_170419_epi23', 'm6268d23_170419_epi26', 'm6233d2_170503_epi22', 'm6233d2_170503_epi25', ...
            'm6233d2_170503_epi28', 'm6233d3_170503_epi21', 'm6233d3_170503_epi26', 'm6233d3_170503_epi29', ...
            'm6233d23_170503_epi20', 'm6233d23_170503_epi24', 'm6233d23_170503_epi31', 'm6233d23_170503_epi33', ...
            'm6233d23_170503_epi34'};
    %}
    %%{
    subjs = {'m5793d3_170418_epi22', 'm5793d3_170418_epi24', 'm5793d23_170418_epi36', 'm6224d23_170411_epi28', 'm6224d3_170411_epi29', ...
             'm6224d2_170411_epi32', 'm6268d2_170419_epi19', 'm6268d2_170419_epi21', 'm6268d23_170419_epi23', 'm6233d23_170503_epi20', ...
             'm6233d3_170503_epi21', 'm6233d3_170503_epi29', 'm6233d23_170503_epi34'};
    %}
    intrp_method = 'nearest';
    load colmap_trP.mat
    startSlice = 7;
    endSlice = 9;
    clstThrld = 4;
    thrldH = 1.7;
    c = 2;
    %%{
    for s = 4:4%length(subjs)
        Vout = spm_vol([path subjs{s} '/anat_img/anatImage.nii']);
        map = spm_read_vols(Vout);
        Vout = spm_vol([path subjs{s} '/epi_data/' maps{c} '.nii']);
        B = spm_read_vols(Vout);
        M = zeros(size(map));
        if size(map, 1) > size(B, 1)
            [x y X Y] = interpAnatFnc(zeros(size(B, 1), size(B, 2)), zeros(size(map, 1), size(map, 2)));
            clstT = clstThrld*((size(map, 1)/size(B, 1))^2);
            for cc = startSlice:endSlice
                P = interp2(x, y, squeeze(B(:, :, cc)), X, Y, intrp_method);
                P(isnan(P)) = 0;
                M(:, :, cc) = P.*voxel_thresholdS(clstT, (P > thrldH));
            end
        end
        if max(M(:)) > 0
            M = M/max(M(:));
            displayFigureR(startSlice, endSlice, map, M, tray, 1);
            title(subjs{s});
        end
        %anat_image = [path subjs{s} '/anat_img/anatImage.nii'];
        %file_name = [path subjs{s} '/epi_data/' maps{c} '.nii'];
        %vizMapH(anat_image, file_name, mp_type, thrldH, thrldL, intp_method, clstThrld, show_slice, subject, contrast, startSlice, endSlice);
    end
    %}
    %{
    subjs = {'m5793', 'm6224', 'm6268', 'm6233'};
    for s = 2:2%length(subjs)
        Vout = spm_vol([path subjs{s} '/anatImage.nii']);
        map = spm_read_vols(Vout);
        Vout = spm_vol([path subjs{s} '/' maps{c} '.nii']);
        B = spm_read_vols(Vout);
        M = zeros(size(map));
        if size(map, 1) > size(B, 1)
            [x y X Y] = interpAnatFnc(zeros(size(B, 1), size(B, 2)), zeros(size(map, 1), size(map, 2)));
            clstT = clstThrld*((size(map, 1)/size(B, 1))^2);
            for cc = startSlice:endSlice
                P = interp2(x, y, squeeze(B(:, :, cc)), X, Y, intrp_method);
                P(isnan(P)) = 0;
                M(:, :, cc) = P.*voxel_thresholdS(clstT, (P > thrldH));
            end
        end
        if max(M(:)) > 0
            M = M/max(M(:));
            displayFigureR(startSlice, endSlice, map, M, tray, 1);
            title(subjs{s});
        end
    end
    %}
    function displayFigureR(startSlice, endSlice, map, M, tray, tp)
        Map = zeros(size(M));
        for slice = startSlice:endSlice
            mP = squeeze(M(:, :, slice));
            d1 = (mP > 0);
            maP = squeeze(map(:, :, slice)).*~d1;
            maP = contrastAdjust(maP);
            maP = 0.48*maP/max(maP(:));
            Map(:, :, slice) = maP + 0.52*d1 + 0.48*mP;
        end
        Map = permute(Map, [2 1 3]);
        disp_Im = [];
        for slice = startSlice:endSlice
            disp_Im = [disp_Im squeeze(Map(:, :, slice))];
        end
        no_slice = endSlice - startSlice + 1;
        if max(disp_Im(:)) > 0.5
            [xsize, ysize] = size(squeeze(Map(:, :, slice)));
            sq_dim = ceil(sqrt(no_slice));
            Im_wdth = sq_dim;
            while (Im_wdth^2 - no_slice >= sq_dim && sq_dim > 1)
                Im_wdth = Im_wdth - 1;
            end
            disIm = zeros(Im_wdth*xsize, sq_dim*ysize);
            figure;
            if Im_wdth ~= 1
                for c2 = 1:Im_wdth - 1
                    disIm((c2 - 1)*xsize + 1:c2*xsize, :) = disp_Im(:, (c2 - 1)*sq_dim*ysize + 1:c2*sq_dim*ysize);
                end
                if sq_dim^2 ~= no_slice
                    disIm(c2*xsize + 1:(c2 + 1)*xsize, 1:ysize*(no_slice - c2*sq_dim)) = disp_Im(:, c2*sq_dim*ysize + 1:end);
                else
                    disIm(c2*xsize + 1:(c2 + 1)*xsize, 1:ysize*sq_dim) = disp_Im(:, c2*sq_dim*ysize + 1:end);
                end
            else
                disIm = disp_Im;
            end
            if tp == 1
                imagesc(disp_Im);
            else
                imagesc(disIm);
            end
            axis off;
            colormap(tray);
        end
    end
    
    %{
    maps = {'beta_0001', 'rspmT_0001', 'PsigC', 'Stim', 'StimH', 'StimT', 'Bstim', 'BstimH', 'BstimT', 'snrRaw', 'Zstim', 'ZstimT'};
    path = '/Users/mishraa/Desktop/arm/limin_obliqLmn/BOLD/HRS_data/';
    subjs = {'m5793d2_170418_epi21', 'm5793d2_170418_epi23', 'm5793d2_170418_epi26', 'm5793d2_170418_epi31', ...
            'm5793d2_170418_epi32', 'm5793d2_170418_epi34', 'm5793d3_170418_epi22', 'm5793d3_170418_epi24', ...
            'm5793d3_170418_epi32', 'm5793d3_170418_epi35', 'm5793d23_170418_epi20', 'm5793d23_170418_epi25', ...
            'm5793d23_170418_epi28', 'm5793d23_170418_epi33', 'm5793d23_170418_epi36', 'm6224d2_170411_epi27', ...
            'm6224d2_170411_epi30', 'm6224d2_170411_epi32', 'm6224d3_170411_epi26', 'm6224d3_170411_epi29', ...
            'm6224d23_170411_epi28', 'm6224d23_170411_epi31', 'm6268d2_170419_epi19', 'm6268d2_170419_epi21', ...
            'm6268d3_170419_epi20', 'm6268d3_170419_epi22', 'm6268d23_170419_epi17', 'm6268d23_170419_epi18', ...
            'm6268d23_170419_epi23', 'm6268d23_170419_epi26', 'm6233d2_170503_epi22', 'm6233d2_170503_epi25', ...
            'm6233d2_170503_epi28', 'm6233d3_170503_epi21', 'm6233d3_170503_epi26', 'm6233d3_170503_epi29', ...
            'm6233d23_170503_epi20', 'm6233d23_170503_epi24', 'm6233d23_170503_epi31', 'm6233d23_170503_epi33', ...
            'm6233d23_170503_epi34'};
    c = 2;
    for s = 1:6%length(subjs)
        Vout = spm_vol([path subjs{s} '/epi_data/' maps{c} '.nii']);
        B = spm_read_vols(Vout);
        Vout = spm_vol([path subjs{1} '/epi_data/' maps{c} '.nii']);
        Vout.fname = [path subjs{s} '/epi_data/' maps{c} '.nii'];
        spm_write_vol(Vout, B);
    end
    %}
        
    
    
    
    
    