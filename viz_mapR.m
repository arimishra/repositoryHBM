%% Program written by Arabinda Mishra, October 20, 2014
clear all; 
close all;
warning off;
    %%{
    subjs = {'m6224_170411_epi16', 'm6224_170411_epi17', 'm6224_170411_epi34', 'm6224_170411_epi35', ...
        'm6224_170411_epi36', 'm6233_170503_epi23', 'm6233_170503_epi27', 'm6233_170503_epi32', ...
        'm6268_170419_epi14', 'm6268_170419_epi15', 'm6268_170419_epi16', 'm6268_170419_epi27', ...
        'm5793_170418_epi17', 'm5793_170418_epi18', 'm5793_170418_epi19'};
    path = '/Users/mishraa/Desktop/arm/limin_obliqLmn/BOLD/HR_data/';
    cvtMap = {'MnCnvt_ar3b_d2T.nii', 'MnCnvt_ar3b_d2B.nii', 'MnCnvt_ar3b_d2L.nii', 'MnCnvt_ar1_T.nii'};
    thrldH = 0.4;
    vxThrld = 2;
    interp_method = 'nearest';
    %% Lr1-4(0.32), Lr2-1(0.22) , Lr3-3(0.18), Ar1Lr1-3(0.27)
    for s = 1:1%length(subjs)
        for c = 1:3
            cnMap = cvtMap{c};
            if cnMap(end - 4) == 'T'
                slice = 9;
            elseif cnMap(end - 4) == 'B'
                slice = 8;
            else
                slice = 7;
            end
            load colmap_trP.mat
            Vout = spm_vol([path subjs{s} '/anat_img/anatImage.nii']);
            A = spm_read_vols(Vout);
            A = squeeze(A(:, :, slice));
            A = 0.5*A/max(A(:));
            Vout = spm_vol([path subjs{s} '/anat_img/' cnMap]);
            T = spm_read_vols(Vout);
            T = squeeze(T(:, :, slice));
            Vout = spm_vol([path subjs{s} '/anat_img/mask.nii']);
            V = spm_read_vols(Vout);
            V = squeeze(V(:, :, slice));
            d1 = (V ~= 0);
            T = d1.*T;
            T(:, 41:end) = 0;
            [x y X Y] = interpAnatFnc(T, A);
            T = interp2(x, y, T, X, Y, interp_method);
            d1 = (T > thrldH);
            [d1] = voxel_thresholdS(vxThrld*round((size(A, 1)/size(V, 1))^2), d1);
            T = d1.*T;
            keeps = find(T > thrldH);
            mxTmap = max(T(:));
            for c = 1:length(keeps)
                T(keeps(c)) = (T(keeps(c)) - thrldH)/(mxTmap -thrldH);
            end
            dsp = A.*~d1 + 0.52*d1 + 0.48*T;
            figure;
            imagesc(permute(dsp(101:400, 101:400), [2, 1]));
            colormap(tray);
            axis off;
        end
    end

        
    
    
    
    
    