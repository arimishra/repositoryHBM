    function cnvt4d_3dICA(pathname, filename)
        Vout = spm_vol([pathname filename]);
        M = spm_read_vols(Vout);
        for i = 1:size(M, 4)
            if i == 1
                V.fname = [pathname 'image' num2str(i, '%03.0f') '.nii'];
                V.mat = Vout(i).mat;
                V.dim = Vout(i).dim;
                V.dt = Vout(i).dt;
                V.pinfo = Vout(i).pinfo;
                V.n = Vout(i).n;
                V.descrip = Vout(i).descrip;
                V.private.dat.fname = V.fname;
                V.private.dat.dim = Vout(i).dim;
                V.private.descrip = Vout(i).descrip;
                V.private.mat = Vout(i).private.mat;
                V.private.mat0 = Vout(i).private.mat0;
                V.private.dat.dtype = Vout(i).private.dat.dtype;
                V.private.dat.offset = Vout(i).private.dat.offset;
                V.private.mat_intent = Vout(i).private.mat_intent;
                V.private.mat0_intent = Vout(i).private.mat0_intent;
                V.private.dat.scl_slope = Vout(i).private.dat.scl_slope;
                V.private.dat.scl_inter = Vout(i).private.dat.scl_inter;
                spm_write_vol(V, squeeze(M(:, :, :, i)));
            else
                V = spm_vol([pathname 'image' num2str(1, '%03.0f') '.nii']);
                V.fname = [pathname 'image' num2str(i, '%03.0f') '.nii'];
                spm_write_vol(V, squeeze(M(:, :, :, i)));
            end
        end
        for c = 1:size(M, 4)
            Vout = spm_vol([pathname 'image' num2str(c, '%03.0f') '.nii']);
            P = spm_read_vols(Vout);
            Vout = spm_vol([pathname 'image001.nii']);
            Vout.fname = [pathname 'image' num2str(c, '%03.0f') '.nii'];
            spm_write_vol(Vout, P);
        end
        disp('Extraction complete');