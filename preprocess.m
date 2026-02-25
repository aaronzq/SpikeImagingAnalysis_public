function [out] = preprocess(roi_name, varargin)

    addpath('./dependencies');
    addpath('./dependencies/EXTRACT-public');
    addpath('./dependencies/EXTRACT-public/EXTRACT');
    addpath('./dependencies/EXTRACT-public/EXTRACT/main_functions');
    addpath('./dependencies/EXTRACT-public/EXTRACT/main_functions/debug_utils');
    addpath('./dependencies/EXTRACT-public/EXTRACT/main_functions/debug_utils/brewermap');
    addpath('./dependencies/EXTRACT-public/EXTRACT/main_functions/helper_functions');
    addpath('./dependencies/EXTRACT-public/EXTRACT/main_functions/solvers');
    addpath('./dependencies/EXTRACT-public/EXTRACT/main_functions/solvers/base solvers for multiple cells');
    addpath('./dependencies/EXTRACT-public/EXTRACT/main_functions/solvers/cell finding solvers');
    addpath('./dependencies/EXTRACT-public/EXTRACT/main_functions/solvers/cell refinement + frr solvers');
    addpath('./dependencies/EXTRACT-public/EXTRACT/modules');
    addpath('./dependencies/EXTRACT-public/EXTRACT/summary functions');
    addpath('./dependencies/EXTRACT-public/External algorithms');
    addpath('./dependencies/EXTRACT-public/External algorithms/Template scripts');
    addpath('./dependencies/EXTRACT-public/External algorithms/parfor_progress');
    addpath('./dependencies/EXTRACT-public/External algorithms/progressbar');
    addpath('./dependencies/EXTRACT-public/External algorithms/simulation_utils');
    addpath('./dependencies/EXTRACT-public/Learning-materials');
    addpath('./dependencies/EXTRACT-public/Learning-materials/Additional demos');
    addpath('./dependencies/EXTRACT-public/Learning-materials/Additional demos/NWB Demos');
    addpath('./dependencies/EXTRACT-public/Learning-materials/Additional demos/Neurofinder_training_0200');
    addpath('./dependencies/EXTRACT-public/Learning-materials/Additional demos/SVD-denoising for low SNR movies');
    addpath('./dependencies/EXTRACT-public/Learning-materials/Additional demos/Simulating movies');
    addpath('./dependencies/EXTRACT-public/Learning-materials/Lecture-Tutorials');
    addpath('./dependencies/EXTRACT-public/Learning-materials/Lecture-Tutorials/Tutorial-1-Starting-code');
    addpath('./dependencies/EXTRACT-public/Learning-materials/Lecture-Tutorials/Tutorial-2-Parallelization');
    addpath('./dependencies/EXTRACT-public/Learning-materials/Lecture-Tutorials/Tutorial-3-Preprocessing');
    addpath('./dependencies/EXTRACT-public/Learning-materials/Lecture-Tutorials/Tutorial-4-Cellfinding');
    addpath('./dependencies/EXTRACT-public/Learning-materials/Lecture-Tutorials/Tutorial-5-Cell-refinement');
    addpath('./dependencies/EXTRACT-public/Learning-materials/Lecture-Tutorials/Tutorial-6-Final-robust-regression');
    addpath('./dependencies/PCAICA');
    addpath('./dependencies/PCAICA/doc');
    addpath('./dependencies/PCAICA/doc/images');
    addpath('./dependencies/altmany-export_fig');
    addpath('./dependencies/altmany-export_fig/.ignore');
    
    addpath('./preprocessing');
    addpath('./preprocessing/demixing');
    addpath('./preprocessing/denoising');
    addpath('./preprocessing/denoising/LOSSDenoising');
    addpath('./preprocessing/denoising/LOSSDenoising/utilities');
    addpath('./preprocessing/denoising/LOSSDenoising/utilities/viewer');
    addpath('./preprocessing/denoising/Spikes_Pipeline_2004');
    addpath('./preprocessing/denoising/Spikes_Pipeline_2004/PROPACK');
    addpath('./preprocessing/denoising/Spikes_Pipeline_2004/utilities');
    addpath('./preprocessing/detrending');
    addpath('./preprocessing/loading');
    addpath('./preprocessing/loading/chunking');
    addpath('./preprocessing/loading/dcimg');
    addpath('./preprocessing/metadata');
    addpath('./preprocessing/motionCorr');
    
    addpath('./analysis');
    addpath('./analysis/cleanCells');
    addpath('./analysis/plots');
    addpath('./analysis/spikeTiming');
    addpath('./analysis/unmixing');
    
    addpath('./utilities');
    addpath('./utilities/paths');
    addpath('./utilities/saving');

    % roi_name = '20260210/obj16X08W_ASAP6c_M1/roi2';
    disp(['Processing /scratch/users/zqwang9/SLab/ASAP6c/' roi_name]);
    h5Path = fullfile('/scratch/users/zqwang9/SLab/ASAP6c', roi_name, 'results/dataset.h5');

    options.run_moco=true;
    options.run_dtr=true;
    options.run_ext=true;
    options.save_h5_range=true;
    options.frame_range=[];

    if nargin>1 % when function received more variable than roi_name
        options=getOptions(options,varargin);
    end

    %%
    load(strrep(h5Path,'dataset.h5','metadata.mat'), 'bpFilter', 'metadata');
    fps = metadata.fps;
    
    if isfile('log.txt')
        delete('log.txt');
    end
    diary log.txt
    parpool('local', str2num(getenv('SLURM_CPUS_PER_TASK')));
    try
        %% Bandpass and motion correction
        if options.run_moco
            bandPassMovieChunk(h5Path, bpFilter);
            path=char(strrep(h5Path,'.h5', '_bp.h5'));
            motionCorr1Movie(path,'nonRigid', false, 'isRawInput', false, 'dcRemoval', false);
        end

        %% Detrend
        if options.run_dtr
            path=char(strrep(h5Path,'.h5', '_bp_moco.h5'));
            if isfile(path)
                detrending(path, 'samplingRate', fps, 'spatialChunk', true);
            else
                disp(['Expecting ' path ' but file is not there.'])
                out = 0;
                return
            end
        end
        %% Extract (demix)
        if options.run_ext       
            path=char(strrep(h5Path,'.h5','_bp_moco_dtr.h5'));
            if isfile(path)
                tic;
                if isempty(options.frame_range)
                    runEXTRACT(path,'polarityGEVI','pos','cellRadius',15,'removeBackground',true,'method','robust');
                else
                    if options.save_h5_range % save motion corrected movie for visual examination
                        h5_target_name = char(strrep(h5Path,'.h5', '_bp_moco.h5'));
                        h5_save_name = char(strrep(h5Path,'.h5', ['_bp_moco' num2str(options.frame_range(1)) '_' num2str(options.frame_range(2)) '.h5']));
                        save_h5_range(h5_target_name, options.frame_range, h5_save_name)               
                    end
                    runEXTRACT(path,'polarityGEVI','pos','cellRadius',15,'removeBackground',true,'method','robust', 'frameRange', options.frame_range);
                end
                toc; 
            else
                disp(['Expecting ' path ' but file is not there.']);
                out = 0;
                return;
            end
        end
        % disp('Preprocess finished.');
        disp(['Finish processing /scratch/users/zqwang9/SLab/ASAP6c/' roi_name]);
        out = 1;

    catch exception
        fprintf('Error occurred: %s\n', exception.message);
        out = 0; % Indicate failure
    end

end