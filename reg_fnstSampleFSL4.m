clear all
close all
warning off
    setenv('FSLDIR', '/usr/local/fsl');
    setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');
    outputPath = '/Users/mishraa/Desktop/arm/limin_obliqLmn/codePackage/sampleData/';
    clear matlabbatch;
    load('prepS1sample.mat');
    spm_jobman('run', matlabbatch);
    im2b_registered = [outputPath 'anat_img/med_img.nii'];
    targetImg = [outputPath 'anat_img/anatImage.nii'];
    regImage = [outputPath 'anat_img/rmed_img.nii.gz'];
    transfMat = [outputPath 'anat_img/trnf.mat'];
    input_Unix = ['/usr/local/fsl/bin/flirt -in ' im2b_registered ' -ref ' targetImg ...
        ' -searchcost normcorr -cost normmi -out ' regImage ' -omat ' transfMat ...
        ' -2D -datatype float -searchrx -5 5 -searchry -5 5 -searchrz 0 0 -interp spline'];
    disp('In Process ....');
    unix(input_Unix);
    gunzip(regImage);
    delete(regImage);
    Vout = spm_vol([outputPath 'prep_rimg4D.nii']);
    M = spm_read_vols(Vout);
    img4_3Dconv([outputPath 'prep_rimg4D.nii']);
    num_vol = size(M, 4);
    for i = 1:num_vol
        im2b_registered = [outputPath 'image' num2str(i, '%03.0f') '.nii'];
        regImage = [outputPath 'rimage' num2str(i, '%03.0f') '.nii.gz'];
        input_Unix = ['/usr/local/fsl/bin/flirt -in ' im2b_registered ' -ref ' targetImg ' -out ' regImage ...
                        ' -2D -applyxfm -init ' transfMat];
        unix(input_Unix);
        gunzip(regImage);
        delete(regImage);
    end
    clear matlabbatch;
    ofuncs = spm_select('FPList', outputPath, '^rimage.*\.nii$');
    matlabbatch{1}.spm.util.cat.name = [outputPath 'prep_rimg4D.nii'];
    matlabbatch{1}.spm.util.cat.vols = cellstr(ofuncs);
    matlabbatch{1}.spm.util.cat.dtype = 4;
    spm_jobman('run', matlabbatch);
    zpF = matlabbatch{1}.spm.util.cat.name;
    delete([zpF(1:end - 4) '.mat']);
    for c = 1:num_vol
        delete([outputPath 'image' num2str(c, '%03.0f') '.nii']);
        delete([outputPath 'rimage' num2str(c, '%03.0f') '.nii']);
    end
    clear matlabbatch;
    load(['prepS2sample.mat']);
    spm_jobman('run', matlabbatch);
        