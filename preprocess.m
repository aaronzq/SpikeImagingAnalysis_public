function [out] = preprocess(h5Path)

    load(strrep(h5Path,'dataset.h5','metadata.mat'), 'bpFilter', 'metadata');
    fps = metadata.fps;
    
    if isfile('log.txt')
        delete('log.txt');
    end
    diary log.txt
    
    try
        %% Bandpass and motion correction
        bandPassMovieChunk(h5Path, bpFilter);
        path=char(strrep(metadata.h5Path,'.h5', '_bp.h5'));
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

end