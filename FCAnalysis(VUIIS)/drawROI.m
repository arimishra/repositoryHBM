%% Program to write mask.
close all;
clear all;
warning off;
clc;
    %% Write ROI mask for any purpose
    %% Written by Arabinda MIshra 01.18.2017
    disp('Select Image: You can either select 3D or 4D EPI image....');
    [filenameF, pathnameF] = uigetfile('*.*', 'Select fMRI data/median Image ... ');
    Vout = spm_vol([pathnameF filenameF]);
    F = spm_read_vols(Vout);
    if numel(size(F)) > 3
        F = permute(median(permute(F, [4 1 2 3])), [2 3 4 1]);
        writeMed_img(pathnameF, filenameF, F);
        Vout = spm_vol([pathnameF 'med_img.nii']);
        M = spm_read_vols(Vout);
    end
    mask = input('Name ROI/Mask File name: ', 's');
    slice_st = str2num(input('Start slice: ', 's'));
    slice_end = str2num(input('Start slice: ', 's'));
    Vm = zeros(size(M));
    VmT = zeros(size(M));
    for slice_no = slice_st:slice_end
        Flag_done = 1;
        Repeat_Flag = 1;
        while (Flag_done == 1)
            dsp = squeeze(M(:, :, slice_no));
            imagesc(flip(permute(dsp, [2 1]), 1));
            colormap(gray);
            axis off;

            if Repeat_Flag == 1
                c1 = 'y';
            end
            if c1 == 'y'
                r = roipoly;
                r = ~(flipdim(permute(r, [2, 1]), 2) & squeeze(M(:, :, slice_no)));
                M1 = r.*squeeze(M(:, :, slice_no)) + ~r*max(M(:));
                close all;
                figure;
                imagesc(flipdim(permute(M1, [2 1]), 1));
                colormap(gray);
                axis off;
                Vm(:, :, slice_no) = ~r*max(M(:));
                Flag_done = 0;
                if slice_no == slice_end
                    Flag_done == 0;
                end
                close all;
            else
                Flag_done = 0;
                close all;
            end
        end
        VmT = VmT + Vm;
    end
    Vout.fname = [pathnameF mask '.nii'];
    spm_write_vol(Vout, VmT);
    close all;
    
    function writeMed_img(pathname, filename, F)
        Vout = spm_vol([pathname filename]);
        V.fname = [pathname 'med_img.nii'];
        V.mat = Vout(1).mat;
        V.dim = Vout(1).dim;
        V.dt = Vout(1).dt;
        V.pinfo = Vout(1).pinfo;
        V.n = Vout(1).n;
        V.descrip = '3D image';
        V.private.dat.fname = V.fname;
        V.private.dat.dim = Vout(1).dim;
        V.private.dat.dtype = Vout(1).private.dat.dtype;
        V.private.dat.offset = Vout(1).private.dat.offset;
        V.private.dat.scl_slope = Vout(1).private.dat.scl_slope;
        V.private.dat.scl_inter = Vout(1).private.dat.scl_inter;
        V.private.mat = Vout(1).private.mat;
        V.private.mat_intent = Vout(1).private.mat_intent;
        V.private.mat0 = Vout(1).private.mat0;
        V.private.mat0_intent = Vout(1).private.mat0_intent;
        V.private.descrip = Vout(1).private.descrip;
        spm_write_vol(V, F);
    end
    
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    
    
    