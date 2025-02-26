%% Steady-state functional connectivity analysis toolbox for SPM12
% the spm_config function that sets up the menus etc in the SPM12 job
% manager.
%
% Authors:
%     Arabinda Mishra
%
% $Id: spm_config_stdcnvty.m 28 2016-04-10 20:43:41Z mishra $
function FncCnvty = tbx_cfg_FncCnvty

    if ~isdeployed 
        addpath(fullfile(spm('dir'), 'toolbox', 'FCAnalysis(VUIIS)'));
    end
    
    %----------------------------------------------------------------------
    roi      = cfg_files;
    roi.name = 'Select ROI image';
    roi.tag  = 'roi';
    roi.filter = 'image';
    roi.num  = [1 100];
    roi.val  = {''};
    roi.help = {'The seed ROI is specified via a binary image mask with the same geometry as ' ...
        'the functional data. This field is mandetory for seed region based connectivity analysis. ' ...
        'The inter ROI correlation matrix is stored as a .mat file in the fmri data folder (intrROIcor.mat) ' ...
        'If you enter for example 10 rois it will be a 9x9 matrix ' ...
        'There are two matrix stored in the folder corMat/_r having inter ROI values calculated differently.'};
    %----------------------------------------------------------------------
    glob      = cfg_menu;
    glob.name = 'Remove global signal?';
    glob.tag  = 'glob';
    glob.labels = {'Yes', 'No'};
    glob.values  = {1, 0};
    glob.val = {0};
    glob.help = {'The whole-volume global signal change may be removed in a regression step, which is ' ...
        'the standard approach following Raichle et al. Not a mandetory fiels however'};

    %----------------------------------------------------------------------
    mot      = cfg_files;
    mot.name = 'Motion parameters?';
    mot.tag  = 'mot';
    mot.filter = '.txt';
    mot.num  = [1 1];
    mot.val = {''};
    mot.help = {'Select the appropriate text file of realignment parameters corresponding to the image data ' ...
        '(rp_motion.txt).  The motion and other confounds will be regressed out from the data prior to connectivity calculation.'};

    %----------------------------------------------------------------------
    conf      = cfg_files;
    conf.name = 'RETROICOR';
    conf.tag  = 'conf';
    conf.filter = '.physiol';
    %conf.filter = '.log';
    %conf.filter = '.mat';
    conf.num = [1 2];
    conf.val = {''};
    conf.help = {'To perform spatial phase correction using physiological signal. Input the physiological parameter file' ...
        'The input signal is contained in file should be saved with extention .txt. Choose Cardiac signal first.'};

    %----------------------------------------------------------------------
    ftsave      = cfg_menu;
    ftsave.name = 'Save filtered files?';
    ftsave.tag  = 'ftsave';
    ftsave.labels = {'Yes', 'No'};
    ftsave.values  = {1, 0};
    ftsave.val  = {0};
    ftsave.help = {'You may save the pre-processed image files with global signal, motion parameters, and additional '...
        'confounds removed, and band pass filter applied for future analysis.'};

    %----------------------------------------------------------------------
    maptyp      = cfg_menu;
    maptyp.name = 'Create Connectivity Map';
    maptyp.tag  = 'maptyp';
    maptyp.labels = {'Mean Connectivity Map', 'Voxelwize Connectivity Map'};
    maptyp.values = {1, 0};
    maptyp.val = {1};
    maptyp.help = {'Connectivity maps may be created based on the mean time series for the seed region, or can be ' ...
        'created for each individual voxel in the ROI.'};

    %----------------------------------------------------------------------
    svbeta      = cfg_menu;
    svbeta.name = 'Save beta value?';
    svbeta.tag  = 'svbeta';
    svbeta.labels = {'Yes', 'No'};
    svbeta.values = {1, 0};
    svbeta.val = {1};
    svbeta.help = {'Functional connectivity is expressed as correlation coefficient, the Fisher-transformed coefficient ' ...
        '(Z-score), and the beta values are an option.'};

    %----------------------------------------------------------------------
    img      = cfg_files;
    img.name = 'Images(.nii/.img)';
    img.tag  = 'img';
    %img.filter = 'image';
    img.filter = '.nii*';
    img.num  = [1 2000];
    img.help = {'Select the preprocessed fMRI time series data for connectivity analysis. '};

    %----------------------------------------------------------------------
    cbr      = cfg_exbranch;
    cbr.name = 'Confounds';
    cbr.tag  = 'cbr';
    cbr.val = {glob, mot, conf};
    cbr.help = {'These signals are removed in a regression step prior to connectivity calculation. ' ...
        'i.e., 1. global signal, 2. Muscle white matter regressor (Can be any confounding parameter) ' ...
        '3. Respiratory and cardiac signal correction using RETROICOR'};

    %----------------------------------------------------------------------
    tr      = cfg_entry;
    tr.name = 'Repetition time (TR)';
    tr.tag = 'tr';
    tr.strtype = 'e';
    tr.num = [1 1];
    tr.val = {0};
    tr.help = {'Enter repetition time in seconds.'};
    
    %----------------------------------------------------------------------
    lpf      = cfg_entry;
    lpf.name = 'BPF cutoff Freq.(0.01 - ?Hz)';
    lpf.tag = 'lpf';
    lpf.strtype = 'e';
    lpf.num = [1 1];
    lpf.val = {0.13};
    lpf.help = {'Default lowpass cutoff frequency is 0.1Hz. Taken into consideration for rsFMRI data analysis.'};
    
    %----------------------------------------------------------------------
    onsetTime = cfg_entry;
    onsetTime.name = 'Onset Time';
    onsetTime.tag = 'onsetTime';
    onsetTime.strtype = 'e';
    onsetTime.num = [1 1];
    onsetTime.val = {0, 0};
    onsetTime.help = {'Onset time for baseline and drug effect signal change calculation.'};
    
    %----------------------------------------------------------------------
    base_window = cfg_entry;
    base_window.name = 'Window width(Baseline/Drug)';
    base_window.tag = 'base_window';
    base_window.strtype = 'e';
    base_window.num = [1 1];
    base_window.val = {0, 0};
    base_window.help = {'Window width (time points) for baseline and drug effect signal change calculation.'};
    
    %----------------------------------------------------------------------
    slide_win      = cfg_menu;
    slide_win.name = 'Need Sliding Window?';
    slide_win.tag  = 'slide_win';
    slide_win.labels = {'Yes', 'No'};
    slide_win.values  = {1, 0};
    slide_win.val = {0};
    slide_win.help = {'Slideing Window saves multiple %Signal change .nii file. Default sliding window moves 5 volumes each ' ...
        '(Continuoos) step. You may choose descrete moves too.'};
    
    %----------------------------------------------------------------------
    move_win      = cfg_menu;
    move_win.name = 'Continuous/Descrete move?';
    move_win.tag  = 'move_win';
    move_win.labels = {'Didcrete', 'Continuous'};
    move_win.values  = {1, 0};
    move_win.val = {0};
    move_win.help = {'Select continuous or descrete window to slide. Default sliding window moves 5 volumes each ' ...
        '(Continuoos) step.'};
    
    %----------------------------------------------------------------------
    pValueTh      = cfg_entry;
    pValueTh.name = 'p Value Threshold';
    pValueTh.tag = 'pValueTh';
    pValueTh.strtype = 'e';
    pValueTh.num = [1 1];
    pValueTh.val = {0.05};
    pValueTh.help = {'Enter Uncorrcted p Value threshold.'};
    
    %----------------------------------------------------------------------
    onset_t      = cfg_entry;
    onset_t.name = 'Onset(block) Time';
    onset_t.tag = 'onset_t';
    onset_t.strtype = 'e';
    onset_t.num = [1 4];
    onset_t.val = {[0 0 0 0]};
    onset_t.help = {'Enter start volume on/off duration and number of volumes eliminated in the beginning of acquisitio.', ...
                    'if on off duration is 30/60s and stimulus occurs after 45s rest at TR = 3, and first 5 volumes in the', ...
                    'begining are exluded enter 5 15 10 20 inplace of zeros. The tool automatically designs the blockand', ...
                    'convolve it with the spm hrf function to locate active regions'};
    
    %----------------------------------------------------------------------
    mask      = cfg_files;
    mask.name = 'External mask';
    mask.tag  = 'mask';
    mask.filter = 'image';
    mask.num  = [1 1];
    mask.val = {''};
    mask.help = {'Optional binary mask image.  If not supplied, an intensity threshold (SPM-antimode) is used for analysis.'};
    
    %----------------------------------------------------------------------
    msk      = cfg_files;
    msk.name = 'External Mask';
    msk.tag  = 'msk';
    msk.filter = 'image';
    msk.num  = [1 1];
    msk.val  = {''};
    msk.help = {'Input ROI .nii files one by one for the Inter ROI connectivity analysis'};
    
    %----------------------------------------------------------------------
    princ      = cfg_files;
    princ.name = 'PCA Analysis(I/p) or WM CSF regression?';
    princ.tag  = 'princ';
    princ.filter = 'image';
    princ.num  = [1 3];
    princ.val = {''};
    princ.help = {'Input binary mask image (.img/.nii) for principal component analysis. Incase WM and CSF need to be ' ...
        'regressed as nuiscance input both masks'};
    
    %----------------------------------------------------------------------
    bld_cbv      = cfg_menu;
    bld_cbv.name = 'BOLD/CBV';
    bld_cbv.tag  = 'bld_cbv';
    bld_cbv.labels = {'BOLD', 'CBV'};
    bld_cbv.values  = {0, 1};
    bld_cbv.val = {0};
    bld_cbv.help = {'Select BOLD/CBV contrast.'};
    
    %----------------------------------------------------------------------
    slide_window      = cfg_entry;
    slide_window.name = 'Sliding Window width(Time Points)';
    slide_window.tag = 'slide_window';
    slide_window.strtype = 'e';
    slide_window.num = [1 1];
    slide_window.val = {0};
    slide_window.help = {'A sliding window connectivity mapping computes the dynamic nature of the connectivity over the whole ' ...
        'brain. User need to specify the width of the window. This is not a mandetory field, if unspecified it considers ' ...
        'no dynamics in the signal.'};
    
    %----------------------------------------------------------------------
    som_maptype        = cfg_menu;
    som_maptype.name   = 'Select Cluster Basis';
    som_maptype.tag    = 'som_maptype';
    som_maptype.labels = {'Connectivity Map', 'Time Course'};
    som_maptype.values = {1, 0};
    som_maptype.val    = {1};
    som_maptype.help   = {'Select connectivity maps as iput to the the clustering analysis. The basis pf the SOM clustering is ' ...
        'the connectivity of each voxel with rest of the brain. However, if you select the set of time series data volumes it ' ...
        'will cluster the time series associated with ROI voxels.'};
    
    %----------------------------------------------------------------------
    dtsom        = cfg_files;
    dtsom.name   = 'Select Voxel-based connectivity Maps(3D/4D';
    dtsom.tag    = 'dtsom';
    dtsom.filter = 'image';
    dtsom.num    = [1 2000];
    dtsom.help   = {'Input the connectivity maps which are supposed to be clustered into different groups',...
    'or you may select the fMRI time series data in order to cluster the time course as the basis of clustering'};

    %----------------------------------------------------------------------
    som_gr         = cfg_entry;
    som_gr.name    = 'Select number of groups';
    som_gr.tag = 'som_gr';
    som_gr.strtype = 'e';
    som_gr.num     = [1 1];
    som_gr.val     = {2};
    som_gr.help    = {'Enter number of groups you like to clissify into:'};
    
    %----------------------------------------------------------------------
    somitr         = cfg_entry;
    somitr.name    = 'Number of iterations(default:240)';
    somitr.tag     = 'somitr';
    somitr.strtype = 'e';
    somitr.num     = [1 1];
    somitr.val     = {240};
    somitr.help    = {'Enter number of iterations SOM has to run: for example 100'};
    
    %----------------------------------------------------------------------
    somsv        = cfg_menu;
    somsv.name   = 'Save connectivity maps?';
    somsv.tag    = 'somsv';
    somsv.labels = {'Yes', 'No'};
    somsv.values = {1, 0};
    somsv.val    = {1};
    somsv.help   = {'Functional connectivity maps clustered in to different groups are saved and the group number can be specified ' ...
        'in the file name By default the clustered maps are not saved in the parent directory.'};
    
    %----------------------------------------------------------------------
    ssc      = cfg_exbranch;
    ssc.name = 'Subject';
    ssc.tag  = 'ssc';
    ssc.val  = {img, roi, onset_t, bld_cbv, slide_window, tr, lpf, cbr, mask, princ, maptyp, ftsave, svbeta, pValueTh};
    ssc.help = {'Steady-state functional connectivity analysis'};
    
    %----------------------------------------------------------------------
    ssdt      = cfg_exbranch;
    ssdt.name = 'Data & Design';
    ssdt.val = {ssc};
    ssdt.help = {'Functional connectivity tool box(VUIIS) designed for resting state and stimulus driven fMRI data analysis ' ...
                'Written implemented on spm12 toolbox by Arabinda Mishra, PhD. '};
    
    %----------------------------------------------------------------------
    ssa      = cfg_exbranch;
    ssa.name = 'Functional Connectivity Analysis';
    ssa.tag  = 'ssa';
    ssa.val  = {ssdt};
    ssa.prog = @tbx_cfg_FncCnvty_func;
    ssa.help = {'Steady state or default mode analysis requires global, spatial motion correction parameters and additional ' ...
        'regressors for preprocessing the data. However, all fields are not mandetory. The user has a choice over selecting ' ...
        'these parameters, you may choose not include any or few of them.'};

    %----------------------------------------------------------------------
    icac      = cfg_exbranch;
    icac.name = 'Subject';
    icac.tag  = 'ssc';
    icac.val  = {img, tr, lpf, onset_t, mask, princ, cbr};
    icac.help = {'Independent Component Analysis'};
    
    %----------------------------------------------------------------------
    icadt        = cfg_exbranch;
    icadt.name   = 'Data & Design';
    icadt.val = {icac};
    icadt.help   = {'Subjects or sessions for which ICA analysis has to be performed.'};
    
    %----------------------------------------------------------------------
    ica      = cfg_exbranch;
    ica.name = 'Pre-Process Data';
    ica.tag  = 'ica';
    ica.val  = {icadt};
    ica.prog = @ICA_analysis;
    ica.help = {'Pre-process data for ICA/Other analysis of data'};
    
    %----------------------------------------------------------------------
    somsbj      = cfg_exbranch;
    somsbj.name = 'Subject';
    somsbj.tag  = 'somsbj';
    somsbj.val  = {img, msk, tr, cbr, princ, mask, slide_window};
    somsbj.help = {'Functional Correlation analysis to follows...'};
    
    %----------------------------------------------------------------------
    somsbj      = cfg_exbranch;
    somsbj.name = 'Subject';
    somsbj.tag  = 'somsbj';
    somsbj.val  = {som_maptype, dtsom, roi, somitr, som_gr, somsv, mask};
    somsbj.help = {'Select connectivity maps for a particular region of interest. Preferably they should be analyzed on voxel basis'};
    
    %----------------------------------------------------------------------
    somdt        = cfg_exbranch;
    somdt.name   = 'Data & Design';
    somdt.val = {somsbj};
    somdt.help   = {'Subjects or sessions for which clustering analysis has to be performed'};
    
    %----------------------------------------------------------------------
    som      = cfg_exbranch;
    som.name = 'SOM Clustering';
    som.tag  = 'som';
    som.val = {somdt};
    som.prog = @som_clust;
    som.help = {'Clustering of subregions in the ROI based on there connectivity maps wrt whole volume'};
    
    %----------------------------------------------------------------------
    FarmCosbj      = cfg_exbranch;
    FarmCosbj.name = 'Subject';
    FarmCosbj.tag  = 'FarmCosbj';
    FarmCosbj.val  = {img, roi, tr, onsetTime, base_window, slide_win, move_win, cbr};
    FarmCosbj.help = {'Select connectivity maps for a particular region of interest. Preferably they should be analyzed on voxel basis'};
    
    %----------------------------------------------------------------------
    FarmCodt        = cfg_exbranch;
    FarmCodt.name   = 'Data & Design';
    FarmCodt.val = {FarmCosbj};
    FarmCodt.help   = {'Subjects or sessions for which Faramacological Analysis to be performed'};
    
    %----------------------------------------------------------------------
    FarmCo      = cfg_exbranch;
    FarmCo.name = 'FarmCo DataAnalysis';
    FarmCo.tag  = 'FarmCo';
    FarmCo.val = {FarmCodt};
    FarmCo.prog = @FarmCo_Analysis;
    FarmCo.help = {'Farmacological Drug response Analysis'};
    
    %----------------------------------------------------------------------
    dynCsbj      = cfg_exbranch;
    dynCsbj.name = 'Subject';
    dynCsbj.tag  = 'ssc';
    dynCsbj.val  = {img, roi};
    dynCsbj.help = {'Select data session'};
    
    %----------------------------------------------------------------------
    dynCdt        = cfg_exbranch;
    dynCdt.name   = 'Data & Design';
    dynCdt.val = {dynCsbj};
    dynCdt.help   = {'Session Data for connectivity dynamics Analysis'};
    
    %----------------------------------------------------------------------
    dynC      = cfg_exbranch;
    dynC.name = 'ConnectivityDynamics(Mazid2011+18)';
    dynC.tag  = 'dynC';
    dynC.val  = {dynCdt};
    dynC.prog = @Dyn_cnvtOptmz;
    dynC.help = {'Functional connectivity dynamics based on Mazid_2011&18'};
    
    %----------------------------------------------------------------------
    FncCnvty        = cfg_choice;
    FncCnvty.name   = 'FCAnalysis(VUIIS)';
    FncCnvty.tag    = 'FCcnvty';
    FncCnvty.values = {ssa, ica};
    %FncCnvty.values = {ica};
    FncCnvty.help   = {'The regression of the whole volume data is performed as a basic preprocessing before the connectivity ' ...
        'analysis. The design matrix contains the data related to motion parameters and artifacts that confounds the fMRI data.' ...
        'The double gamma fitted homodynamic responce function are convolved with the block design data and are modulated with ' ...
        'the seed time course. The connectivity map are the coefficient of regression when the whole volume is analysed wrt. the ' ...
        'seed time course. The final maps can be stored in two different ways, i.e., i)connectivity map wrt mean signal in ' ...
        'the ROI or ii) individual maps'};
    %----------------------------------------------------------------------
    function tbx_cfg_FncCnvty_func(job)
        ssc = job.generic.ssc;
        spm_FCCnvty(ssc.img, ssc.roi, ssc.tr, ssc.lpf, ssc.cbr, ssc.maptyp, ssc.ftsave, ssc.svbeta, ssc.mask, ssc.princ, ...
                    ssc.onset_t, ssc.bld_cbv, ssc.slide_window, ssc.pValueTh);
        
    function ICA_analysis(job)
        ssc = job.generic.ssc;
        spm_ICAPrep(ssc.img, ssc.tr, ssc.lpf, ssc.cbr, ssc.mask, ssc.princ, ssc.onset_t);
    

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
 
    
    
    
    
