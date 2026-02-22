% function [out] = preprocess(h5Path)
    addpath('SpikeImagingAnalysis_public');
    addpath('SpikeImagingAnalysis_public/dependencies');
    addpath('SpikeImagingAnalysis_public/dependencies/PCAICA');
    addpath('SpikeImagingAnalysis_public/dependencies/PCAICA/doc');
    addpath('SpikeImagingAnalysis_public/dependencies/PCAICA/doc/images');
    addpath('SpikeImagingAnalysis_public/dependencies/altmany-export_fig');
    addpath('SpikeImagingAnalysis_public/dependencies/altmany-export_fig/.ignore');
    addpath('SpikeImagingAnalysis_public/preprocessing');
    addpath('SpikeImagingAnalysis_public/preprocessing/demixing');
    addpath('SpikeImagingAnalysis_public/preprocessing/denoising');
    addpath('SpikeImagingAnalysis_public/preprocessing/denoising/LOSSDenoising');
    addpath('SpikeImagingAnalysis_public/preprocessing/denoising/LOSSDenoising/utilities');
    addpath('SpikeImagingAnalysis_public/preprocessing/denoising/LOSSDenoising/utilities/viewer');
    addpath('SpikeImagingAnalysis_public/preprocessing/denoising/Spikes_Pipeline_2004');
    addpath('SpikeImagingAnalysis_public/preprocessing/denoising/Spikes_Pipeline_2004/PROPACK');
    addpath('SpikeImagingAnalysis_public/preprocessing/denoising/Spikes_Pipeline_2004/utilities');
    addpath('SpikeImagingAnalysis_public/preprocessing/detrending');
    addpath('SpikeImagingAnalysis_public/preprocessing/loading');
    addpath('SpikeImagingAnalysis_public/preprocessing/loading/chunking');
    addpath('SpikeImagingAnalysis_public/preprocessing/loading/dcimg');
    addpath('SpikeImagingAnalysis_public/preprocessing/metadata');
    addpath('SpikeImagingAnalysis_public/preprocessing/motionCorr');
    addpath('SpikeImagingAnalysis_public/analysis');
    addpath('SpikeImagingAnalysis_public/analysis/cleanCells');
    addpath('SpikeImagingAnalysis_public/analysis/plots');
    addpath('SpikeImagingAnalysis_public/analysis/spikeTiming');
    addpath('SpikeImagingAnalysis_public/analysis/unmixing');
    addpath('SpikeImagingAnalysis_public/utilities');
    addpath('SpikeImagingAnalysis_public/utilities/paths');
    addpath('SpikeImagingAnalysis_public/utilities/saving');

    roiName = '20260210/obj16X08W_ASAP6c_M1/roi2';
    h5Path = fullfile('/scratch/users/zqwang9/SLab/ASAP6c', roiName, 'results/dataset.h5');

    load(strrep(h5Path,'dataset.h5','metadata.mat'), 'bpFilter', 'metadata');
    fps = metadata.fps;
    
    if isfile('log.txt')
        delete('log.txt');
    end
    diary log.txt
    parpool('local', str2num(getenv('SLURM_CPUS_PER_TASK')));
    try
        %% Bandpass and motion correction
        bandPassMovieChunk(h5Path, bpFilter);
        path=char(strrep(h5Path,'.h5', '_bp.h5'));
        motionCorr1Movie(path,'nonRigid', false,'isRawInput',false,'dcRemoval',false);
    
        %% Detrend
        path=char(strrep(h5Path,'.h5', '_bp_moco.h5'));
        detrending(path, 'samplingRate', fps,'spatialChunk', true);
    
        %% Extract (demix)
        path=char(strrep(h5Path,'.h5','_bp_moco_dtr.h5'));
        tic;runEXTRACT(path,'polarityGEVI','pos','cellRadius',15,'removeBackground',true,'method','robust');toc;  

        disp('Preprocess finished.');
        out = 1;

    catch exception
        fprintf('Error occurred: %s\n', exception.message);
        out = 0; % Indicate failure
    end

% end