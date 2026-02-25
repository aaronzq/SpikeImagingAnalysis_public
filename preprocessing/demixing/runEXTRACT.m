function [output]=runEXTRACT(h5Path,varargin)


%% OPTIONS

options.diary=true;
options.plotFigure=true;
options.diary_name='diary';
options.verbose=true;
options.frameRange=[];
options.polarityGEVI='pos'; % 'neg' 'dual'
options.rank=100;
options.binning=[];
options.partition=[1 1];
options.cellRadius=20;
options.removeBackground=true;
options.method='L2'; %'robust'

%% UPDATE OPTIONS
if nargin>1
    options=getOptions(options,varargin);
end

%% LOAD RAW DCMIG DATA
disps('Starting Neuron Demixing (EXTRACT)')


if ischar(h5Path)
    [folderPath,fname,ext]=fileparts(h5Path);
    if strcmpi(ext,'.h5')||strcmpi(ext,'.hdf5')
        disps('h5 file detected')
    else
        error('not a h5 or Tiff file')
    end
    
    meta=h5info(h5Path);
    dim=meta.Datasets.Dataspace.Size;
    mx=dim(1);my=dim(2);numFrame=dim(3);
    if isempty(options.frameRange)
        options.frameRange=[1 numFrame];
    end
    dataset=strcat(meta.Name,meta.Datasets.Name);
    outputdir=fullfile(folderPath,'DemixingEXTRACT',fname);
    M=h5read(h5Path,dataset,[1 1 options.frameRange(1)],[mx my diff(options.frameRange)+1]);
    disp('Movie succefully loaded')
    if ~isempty(options.binning)
        disps('binning the movie')
        M=imresize3(M,[mx/options.binning my/options.binning diff(options.frameRange)+1],'box');
    end
    if ~exist(outputdir,'dir')
        mkdir(outputdir)
    end
    
else
    M=h5Path;
end
%%
%Initialize config
config=[];
config = get_defaults(config);

%Set some important settings
config.use_gpu=1;
config.verbose = 2;
config.dendrite_aware=0; %If you want dendrites from the movie
config.use_gpu=1;
% Essentials, without these EXTRACT will give an error:
config.avg_cell_radius=options.cellRadius; %Pick a reasonable cell radius

%Optionals, but strongly advised to handpick:
%Movie is small enough that EXTRACT will not automatically partition this,
%but still a good idea to keep these in sight!
config.num_partitions_x=options.partition(1);
config.num_partitions_y=options.partition(2);
% config.temporal_denoising=0;
% config.spatial_highpass_spatial=inf;
% config.cellfind_filter_type='butter'; % Feel free to use your own filters,

% set algo iteration for finding cells
config.cellfind_max_steps=100;

% set algo iteration for cleaning duplicate cells
config.max_iter=5;

% Voltage specific configs
config.trace_output_option='no_constraint';
config.cellfind_min_snr=0;
config.thresholds.T_min_snr=0;
config.cellfind_kappa_std_ratio=0.7; % for final cell find

% my movie are already preprocessed (moco and detrending)
config.preprocess=1;

% parameter to use for Extract paper
if strcmpi(options.method,'L2')   
    config.kappa_std_ratio=100; % L2 solver - used in Science paper    
elseif strcmpi(options.method,'robust')    
    config.kappa_std_ratio=0.7; % Extract solver % SH-20230116
else
    error('requested demixing method not found...');
end

config.background_removal=options.removeBackground;

%Perform the extraction
switch options.polarityGEVI
    case 'pos'
        
        % M=M-1; % -1 added to compute (f-f0)/f0 instead of f/f0
        output=extractor(M,config);
        
        if ischar(h5Path)
            if isempty(options.frameRange)
                fileName=['extractOutput_' options.polarityGEVI '_all' '_' datestr(datetime('now'),'yyyymmddTHHMMSS') '.mat'];
            else
                fileName=['extractOutput_' options.polarityGEVI '_' num2str(options.frameRange(1)) '_' num2str(options.frameRange(2)) '_' datestr(datetime('now'),'yyyymmddTHHMMSS') '.mat'];
            end
            savePath=fullfile(outputdir,fileName);
            save(savePath,'output');
        end
        
    case 'neg'
        
        disps('matrix inversion')
        M=-M+2*mean(M,3); % -1 added to compute (f-f0)/f0 instead of f/f0
%         plot(getPointProjection(M))
        output=extractor(M,config);
        
        if ischar(h5Path)
            fileName=['extractOutput_' options.polarityGEVI '_' options.method '_' datestr(datetime('now'),'yyyymmddTHHMMSS') '.mat'];
            savePath=fullfile(outputdir,fileName);
            save(savePath,'output');
        end
        
    case 'dual'
        
        disps('running Extract twice for pos&neg polarities')
        disps('running positive Extract')
        
        output=extractor(M,config);
        
        options.polarityGEVI='pos';
        if ischar(h5Path)
            fileName=['extractOutput_' options.polarityGEVI '_' options.method '_' datestr(datetime('now'),'yyyymmddTHHMMSS') '.mat'];
            savePath=fullfile(outputdir,fileName);
            save(savePath,'output');
        end
        
        disps('running negative Extract')        
        disps('matrix inversion')
        M=-M+2*mean(M,3); % -1 added to compute (f-f0)/f0 instead of f/f0
        
        output=extractor(M,config);
        
        options.polarityGEVI='neg';
        if ischar(h5Path)
            fileName=['extractOutput_' options.polarityGEVI '_' options.method '_' datestr(datetime('now'),'yyyymmddTHHMMSS') '.mat'];
            savePath=fullfile(outputdir,fileName);
            save(savePath,'output');
        end
end

    function disps(string) %overloading disp for this function
        if options.verbose
            fprintf('%s Demixing EXTRACT: %s\n', datetime('now'),string);
        end
    end
end