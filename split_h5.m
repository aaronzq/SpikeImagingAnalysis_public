%% Separate h5 dataset in n parts for imagej visualization

mainPath = 'C:\Users\Z\Documents\SLab\20260210\obj16X08W_ASAP6c_M1\roi2\results';
h5Name = 'dataset_bp_moco.h5';
savePath = fullfile(mainPath, 'visualization');

if ~exist(savePath, "dir")
    mkdir(savePath);
end

nParts = 2;
chunkSize=1000;
%%
movieDims = h5info(fullfile(mainPath, h5Name)).Datasets.Dataspace.Size;
height = movieDims(1);
width = movieDims(2);
num_pixels = height * width;
numFrame = movieDims(3);
numFramePart = floor(numFrame/nParts);
partFirstLast = chunkFrames(numFramePart,[1 numFrame]);

get_frames = @(x) h5read(fullfile(mainPath, h5Name), '/mov', [1 1 x(1)], [movieDims(1) movieDims(2) x(end)-x(1)+1]);
disp('Splitting big h5 file into parts...')
for i = 1:size(partFirstLast,1)
    h5SaveName = char(strrep(h5Name,'.h5',['_' num2str(i) '.h5']));
    
    chunksFirstLast = chunkFrames(chunkSize,[partFirstLast(i,1) partFirstLast(i,2)]);


    if isfile(fullfile(savePath, h5SaveName))
        disps(['Already found h5 file:' h5SaveName 'deleting!']);
        delete(fullfile(savePath, h5SaveName));
    end
    
    for ichunk=1:size(chunksFirstLast,1) 
        movieBatch = get_frames([chunksFirstLast(ichunk,1),chunksFirstLast(ichunk,2)]);
        
        h5append(fullfile(savePath, h5SaveName), single(movieBatch)); 
        
    end

end
disp('Split done.')