%% GLM for small volume correction

% Subject info.
SN = {};
SN{end+1} = '01'; % PSY
SN{end+1} = '03'; % KIS
SN{end+1} = '04'; % JBH
SN{end+1} = '05'; % HJH %%% check alignment for sub4,5 - 6 deg mask 181210
SN{end+1} = '06'; % KB
SN{end+1} = '07'; % SHY
SN{end+1} = '08'; % LSY
SN{end+1} = '09'; % HJH2
SN{end+1} = '10'; % CES
SN{end+1} = '11'; % YYH
SN{end+1} = '12'; % LHB
SN{end+1} = '13'; % JHY
SN{end+1} = '14'; % KYJ

%% path 
root_dir = '/group_hpc/WMShimLab2/PSY_Color/Colorv3';

%% initialize spm defaults
spm('defaults', 'FMRI');
    spm_jobman('initcfg');

     %% text files
     color_onset = importdata('/group_hpc/WMShimLab/PSY_ColorStudy/loc_color.txt');
     fix_onset = importdata('/group_hpc/WMShimLab/PSY_ColorStudy/loc_fix.txt');
    
for xSN = 3:length(SN)
clear batch;
    
    loc_file_dir = fullfile(root_dir, SN{xSN}, 'Img_data', sprintf('%s_loc.results', SN{xSN}));
    full_mask = fullfile(loc_file_dir, sprintf('full_mask.%s_loc.nii, 1', SN{xSN}));
    glm_output_dir = fullfile(root_dir, SN{xSN}, 'Img_data', sprintf('%s_spm_loc.results', SN{xSN}));
    if ~exist(glm_output_dir,'dir'); mkdir(glm_output_dir); end
    
    %-----------------------------------------------------------------------
    % Job saved on 06-May-2020 19:26:06 by cfg_util (rev $Rev: 6460 $)
    % spm SPM - SPM12 (6906)
    % cfg_basicio BasicIO - Unknown
    %-----------------------------------------------------------------------
    batch{1}.spm.stats.fmri_spec.dir = {glm_output_dir};
    batch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    batch{1}.spm.stats.fmri_spec.timing.RT = 2;
    batch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    batch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
    
    
    %% RUN 1
    ff =[]; xRun = 1;
    loc_run = sprintf('pb01.%s_loc.r%02d.sc_dt_hp_am_blur.nii', SN{xSN}, xRun);
    ff = [ff; spm_select('ExtFPList', loc_file_dir, loc_run, Inf)];
    batch{1}.spm.stats.fmri_spec.sess(1).scans = cellstr(ff);
    %%
    batch{1}.spm.stats.fmri_spec.sess(1).cond(1).name = 'color';
    batch{1}.spm.stats.fmri_spec.sess(1).cond(1).onset = color_onset(xRun,:);
    batch{1}.spm.stats.fmri_spec.sess(1).cond(1).duration = 12;
    batch{1}.spm.stats.fmri_spec.sess(1).cond(1).tmod = 0;
    batch{1}.spm.stats.fmri_spec.sess(1).cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
    batch{1}.spm.stats.fmri_spec.sess(1).cond(1).orth = 1;
    batch{1}.spm.stats.fmri_spec.sess(1).cond(2).name = 'fix';
    %%
    batch{1}.spm.stats.fmri_spec.sess(1).cond(2).onset = fix_onset(xRun,:);
    %%
    batch{1}.spm.stats.fmri_spec.sess(1).cond(2).duration = 12;
    batch{1}.spm.stats.fmri_spec.sess(1).cond(2).tmod = 0;
    batch{1}.spm.stats.fmri_spec.sess(1).cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
    batch{1}.spm.stats.fmri_spec.sess(1).cond(2).orth = 1;
    batch{1}.spm.stats.fmri_spec.sess(1).multi = {''};
    batch{1}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {});
    batch{1}.spm.stats.fmri_spec.sess(1).multi_reg = {''};
    batch{1}.spm.stats.fmri_spec.sess(1).hpf = 128;
    
    %% RUN 2
    ff =[]; xRun = 2;
    loc_run = sprintf('pb01.%s_loc.r%02d.sc_dt_hp_am_blur.nii', SN{xSN}, xRun);
    ff = [ff; spm_select('ExtFPList', loc_file_dir, loc_run, Inf)];
    batch{1}.spm.stats.fmri_spec.sess(2).scans = cellstr(ff);
    
    %%
    batch{1}.spm.stats.fmri_spec.sess(2).cond(1).name = 'color';
    batch{1}.spm.stats.fmri_spec.sess(2).cond(1).onset = color_onset(xRun,:);
    batch{1}.spm.stats.fmri_spec.sess(2).cond(1).duration = 12;
    batch{1}.spm.stats.fmri_spec.sess(2).cond(1).tmod = 0;
    batch{1}.spm.stats.fmri_spec.sess(2).cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
    batch{1}.spm.stats.fmri_spec.sess(2).cond(1).orth = 1;
    batch{1}.spm.stats.fmri_spec.sess(2).cond(2).name = 'fix';
    %%
    batch{1}.spm.stats.fmri_spec.sess(2).cond(2).onset = fix_onset(xRun,:);
    %%
    batch{1}.spm.stats.fmri_spec.sess(2).cond(2).duration = 12;
    batch{1}.spm.stats.fmri_spec.sess(2).cond(2).tmod = 0;
    batch{1}.spm.stats.fmri_spec.sess(2).cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
    batch{1}.spm.stats.fmri_spec.sess(2).cond(2).orth = 1;
    batch{1}.spm.stats.fmri_spec.sess(2).multi = {''};
    batch{1}.spm.stats.fmri_spec.sess(2).regress = struct('name', {}, 'val', {});
    batch{1}.spm.stats.fmri_spec.sess(2).multi_reg = {''};
    batch{1}.spm.stats.fmri_spec.sess(2).hpf = 128;
    
    %% RUN 3
    ff =[]; xRun = 3;
    loc_run = sprintf('pb01.%s_loc.r%02d.sc_dt_hp_am_blur.nii', SN{xSN}, xRun);
    ff = [ff; spm_select('ExtFPList', loc_file_dir, loc_run, Inf)];
    batch{1}.spm.stats.fmri_spec.sess(3).scans = cellstr(ff);
    
    %%
    batch{1}.spm.stats.fmri_spec.sess(3).cond(1).name = 'color';
    batch{1}.spm.stats.fmri_spec.sess(3).cond(1).onset = color_onset(xRun,:);
    batch{1}.spm.stats.fmri_spec.sess(3).cond(1).duration = 12;
    batch{1}.spm.stats.fmri_spec.sess(3).cond(1).tmod = 0;
    batch{1}.spm.stats.fmri_spec.sess(3).cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
    batch{1}.spm.stats.fmri_spec.sess(3).cond(1).orth = 1;
    batch{1}.spm.stats.fmri_spec.sess(3).cond(2).name = 'fix';
    %%
    batch{1}.spm.stats.fmri_spec.sess(3).cond(2).onset = fix_onset(xRun,:);
    %%
    batch{1}.spm.stats.fmri_spec.sess(3).cond(2).duration = 12;
    batch{1}.spm.stats.fmri_spec.sess(3).cond(2).tmod = 0;
    batch{1}.spm.stats.fmri_spec.sess(3).cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
    batch{1}.spm.stats.fmri_spec.sess(3).cond(2).orth = 1;
    batch{1}.spm.stats.fmri_spec.sess(3).multi = {''};
    batch{1}.spm.stats.fmri_spec.sess(3).regress = struct('name', {}, 'val', {});
    batch{1}.spm.stats.fmri_spec.sess(3).multi_reg = {''};
    batch{1}.spm.stats.fmri_spec.sess(3).hpf = 128;
    
    %% RUN 4
    ff =[]; xRun = 4;
    loc_run = sprintf('pb01.%s_loc.r%02d.sc_dt_hp_am_blur.nii', SN{xSN}, xRun);
    ff = [ff; spm_select('ExtFPList', loc_file_dir, loc_run, Inf)];
    batch{1}.spm.stats.fmri_spec.sess(4).scans = cellstr(ff);
    
    %%
    batch{1}.spm.stats.fmri_spec.sess(4).cond(1).name = 'color';
    batch{1}.spm.stats.fmri_spec.sess(4).cond(1).onset = color_onset(xRun,:);
    batch{1}.spm.stats.fmri_spec.sess(4).cond(1).duration = 12;
    batch{1}.spm.stats.fmri_spec.sess(4).cond(1).tmod = 0;
    batch{1}.spm.stats.fmri_spec.sess(4).cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
    batch{1}.spm.stats.fmri_spec.sess(4).cond(1).orth = 1;
    batch{1}.spm.stats.fmri_spec.sess(4).cond(2).name = 'fix';
    %%
    batch{1}.spm.stats.fmri_spec.sess(4).cond(2).onset = fix_onset(xRun,:);
    %%
    batch{1}.spm.stats.fmri_spec.sess(4).cond(2).duration = 12;
    batch{1}.spm.stats.fmri_spec.sess(4).cond(2).tmod = 0;
    batch{1}.spm.stats.fmri_spec.sess(4).cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
    batch{1}.spm.stats.fmri_spec.sess(4).cond(2).orth = 1;
    batch{1}.spm.stats.fmri_spec.sess(4).multi = {''};
    batch{1}.spm.stats.fmri_spec.sess(4).regress = struct('name', {}, 'val', {});
    batch{1}.spm.stats.fmri_spec.sess(4).multi_reg = {''};
    batch{1}.spm.stats.fmri_spec.sess(4).hpf = 128;
    batch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    batch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    batch{1}.spm.stats.fmri_spec.volt = 1;
    batch{1}.spm.stats.fmri_spec.global = 'None';
    batch{1}.spm.stats.fmri_spec.mthresh = 0.8;
    batch{1}.spm.stats.fmri_spec.mask = cellstr(full_mask);
    batch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    batch{1}.spm.stats.review.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    batch{1}.spm.stats.review.display.matrix = 1;
    batch{1}.spm.stats.review.print = 'ps';
    batch{1}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    batch{1}.spm.stats.fmri_est.write_residuals = 1;
    batch{1}.spm.stats.fmri_est.method.Classical = 1;
    batch{1}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    batch{1}.spm.stats.con.consess{1}.tcon.name = 'color>fix';
    batch{1}.spm.stats.con.consess{1}.tcon.weights = [0.25 -0.25 0.25 -0.25 0.25 -0.25 0.25 -0.25];
    batch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    batch{1}.spm.stats.con.delete = 0;
    
    
end

