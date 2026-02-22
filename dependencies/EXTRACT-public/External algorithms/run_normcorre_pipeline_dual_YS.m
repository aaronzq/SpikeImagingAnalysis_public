function run_normcorre_pipeline_dual_YS(input,output,input2,output2,config)
% any additional params can be set in config

nt_template     = config.nt_template;
template        = config.template;
numFrame        = config.numFrame;
nonrigid_mc     = config.nonrigid_mc;
ns_nonrigid     = config.ns_nonrigid;
bandpass        = config.bandpass;
avg_cell_radius = config.avg_cell_radius;
use_gpu         = config.use_gpu;
file_type       = config.file_type;
mask            = config.mask;
mc_template     = config.mc_template;
get_mask        = config.get_mask;
if isfield(config,'grid_size') && ns_nonrigid==1
    nonrigid_mc = config.grid_size 
    disp('      -----CAUTION! Non-rigid motion correction ongoing-----')
else
    disp('      -----CAUTION! Rigid motion correction ongoing-----')
end
if isfield(config,'max_shift')
    max_shift = config.max_shift 
else
    max_shift = round(ns_nonrigid/5);
end

switch file_type
    case 'h5'
        [input_filename,input_datasetname] = parse_movie_name(input);
        movie_info = h5info(input_filename,input_datasetname);
        movie_size = num2cell(movie_info.Dataspace.Size);
        [nx, ny, totalnum] = deal(movie_size{:});
        %
        [input_filename2,input_datasetname2] = parse_movie_name(input2);
        movie_info2 = h5info(input_filename2,input_datasetname2);
        movie_size2 = num2cell(movie_info2.Dataspace.Size);
        [nx2, ny2, totalnum2] = deal(movie_size2{:});                
    case 'tif'
        tiff_info = imfinfo(input);
        nx = tiff_info(1).Height;
        ny = tiff_info(1).Width;
        totalnum = size(tiff_info, 1);
        %
        tiff_info2 = imfinfo(input2);
        nx2 = tiff_info2(1).Height;
        ny2 = tiff_info2(1).Width;
        totalnum2 = size(tiff_info2, 1);        
    case 'tiff'
        tiff_info = imfinfo(input);
        nx = tiff_info(1).Height;
        ny = tiff_info(1).Width;
        totalnum = size(tiff_info, 1);
        %
        tiff_info2 = imfinfo(input2);
        nx2 = tiff_info2(1).Height;
        ny2 = tiff_info2(1).Width;
        totalnum2 = size(tiff_info2, 1);
end

[output_filename,output_datasetname] = parse_movie_name(output);
[output_filename2,output_datasetname2] = parse_movie_name(output2);

if isfile(output_filename)
    delete(output_filename);
end
if isfile(output_filename2)
    delete(output_filename2);
end

count = 0;
chunk_time_size = numFrame;
while count == 0
    try
        h5create(output_filename,output_datasetname,[nx ny totalnum],'Datatype','single','ChunkSize',[nx,ny,chunk_time_size]);
        h5create(output_filename2,output_datasetname2,[nx2 ny2 totalnum2],'Datatype','single','ChunkSize',[nx2,ny2,chunk_time_size]);
        count = 1;
    catch 
        chunk_time_size = round(chunk_time_size/2);
    end
end

windowsize = min(totalnum, numFrame);

startno = [1:windowsize:totalnum];

if numel(startno)>1
    % handling the irregular framenumbers 
    perframes = ones(numel(startno),1)*numFrame;
    lastframes = mod(totalnum,numFrame);
    if lastframes > 0
        perframes(end-1) = perframes(end-1) + lastframes;
        startno(end) = [];
    end
else
    perframes = totalnum;
end

if nt_template > totalnum
    nt_template = totalnum;
end

if isempty(template)
    % create the template FOV
    disp(sprintf('%s: Getting the template for motion correction from the first %s frames', datestr(now), num2str(nt_template) ))
    
    switch file_type
        case 'h5'
            im1 = single(h5read(input_filename, input_datasetname, [1, 1, 1], [nx, ny, nt_template]));   
            im12 = single(h5read(input_filename2, input_datasetname2, [1, 1, 1], [nx2, ny2, nt_template]));   
        case 'tif'
            im1 = single(read_from_tif(input,1,nt_template));
            im12 = single(read_from_tif(input2,1,nt_template));
        case 'tiff'
            im1 = single(read_from_tif(input,1,nt_template));
            im12 = single(read_from_tif(input2,1,nt_template));
    end

    if bandpass
        im1 = spatial_bandpass(im1,avg_cell_radius,10,2,use_gpu);
        im12 = spatial_bandpass(im12,avg_cell_radius,10,2,use_gpu);
    end

    im1_nonrig = im1;
    im1_nonrig2 = im12;

    if get_mask
        figure
        imshow(mean(im1,3),[])
        rect = getrect
        y_min = round(rect(1));
        x_min = round(rect(2));
        y_max = round(rect(1)+rect(3));
        x_max = round(rect(2)+rect(4));
        mask = zeros(nx,ny);
        mask(x_min:x_max,y_min:y_max) = 1;

        figure
        imshow(mean(im12,3),[])
        rect2 = getrect
        y_min2 = round(rect2(1));
        x_min2 = round(rect2(2));
        y_max2 = round(rect2(1)+rect2(3));
        x_max2 = round(rect2(2)+rect2(4));
        mask2 = zeros(nx2,ny2);
        mask2(x_min2:x_max2,y_min2:y_max2) = 1;    
    end
    
    if ~isempty(mask)
        idx_2 = find(sum(mask,1)==0);
        idx_1 = find(sum(mask,2)==0);
        im1(idx_1(:),:,:)=[];
        im1(:,idx_2(:),:)=[];
        
        im12(idx_1(:),:,:)=[];
        im12(:,idx_2(:),:)=[];        
    end

    if mc_template
        options_template = NoRMCorreSetParms('d1',size(im1,1),'d2',size(im1,2),'max_shift',30,'us_fac',30,'grid_size',[size(im1,1),size(im1,2)],'print_msg',0); 
        % The grid size is full FOV
        [~,shifts_template,~,options_tempate] = normcorre_batch(im1,options_template,mean(im1(:,:,1:max(round(nt_template/10),2)),3)) ;
        im1 = apply_shifts(im1,shifts_template,options_template);
        im12 = apply_shifts(im12,shifts_template,options_template);

        if ~isempty(mask)
            options_template.d1=size(im1_nonrig,1);
            options_template.d2=size(im1_nonrig,2);
            options_template.grid_size = [size(im1_nonrig,1),size(im1_nonrig,2),1];
            im1_nonrig = apply_shifts(im1_nonrig,shifts_template,options_template);
            im1_nonrig2 = apply_shifts(im1_nonrig2,shifts_template,options_template);
        else
            im1_nonrig = im1;
            im1_nonrig2 = im12;
        end
        clear options_template
        clear shifts_template
    end

    %im_ds= single(max(im1,[],3));
    im_ds= single(mean(im1,3));
    im_ds_nonrig = single(mean(im1_nonrig,3));
    im_ds2= single(mean(im12,3));
    im_ds_nonrig2 = single(mean(im1_nonrig2,3));    
    figure
    imshow(im_ds,[]);
    drawnow;
    clear im1 im12

else
    disp(sprintf('%s: Template image already assigned.', datestr(now)))
    im_ds = single(template);
    im_ds_nonrig = im_ds;
    
    if ~isempty(mask)
        idx_2 = find(sum(mask,1)==0);
        idx_1 = find(sum(mask,2)==0);
        im_ds(idx_1(:),:)=[];
        im_ds(:,idx_2(:))=[];
    end
end

% motion correction starts below
disp(sprintf('%s: Running motion correction split into %s movies', datestr(now),num2str(numel(startno)) ))

for i=1:numel(startno)
    % run motion correction in each movie block
    disp(sprintf('\t %s: Processing %s/%s movies, %s', datestr(now),num2str(i),num2str(numel(startno)),cd ))

    switch file_type
        case 'h5'
            M = single(h5read(input_filename,input_datasetname,[1,1,startno(i)],[nx,ny,perframes(i)]));
            M2 = single(h5read(input_filename2,input_datasetname2,[1,1,startno(i)],[nx2,ny2,perframes(i)]));
        case 'tif'
            M = single(read_from_tif(input,startno(i),perframes(i)));
            M2 = single(read_from_tif(input2,startno(i),perframes(i)));
        case 'tiff'
            M = single(read_from_tif(input,startno(i),perframes(i)));
            M2 = single(read_from_tif(input2,startno(i),perframes(i)));
    end
    
    M_proc = M;
    M_proc2 = M2;
    if bandpass
        M_proc = spatial_bandpass(M_proc,avg_cell_radius,10,2,use_gpu);
        M_proc2 = spatial_bandpass(M_proc2,avg_cell_radius,10,2,use_gpu);
    end
    if ~isempty(mask)
        M_proc(idx_1(:),:,:)=[];
        M_proc(:,idx_2(:),:)=[];
        M_proc2(idx_1(:),:,:)=[];
        M_proc2(:,idx_2(:),:)=[];        
    end

    % rigid motion correction starts below -------------------------
    disp(sprintf('\t \t %s: Starting rigid motion correction ', datestr(now)))

    options_rigid = NoRMCorreSetParms('d1',size(M_proc,1),'d2',size(M_proc,2),'max_shift',max_shift,'us_fac',30,'grid_size',[size(M_proc,1),size(M_proc,2)],'print_msg',0); 
    % The grid size is full FOV
    [~,shifts_rigid,~,options_rigid] = normcorre_batch(M_proc,options_rigid,im_ds);

    clear M_proc M_proc2
    
    if ~isempty(mask)
        options_rigid.d1=nx;
        options_rigid.d2=ny;
        options_rigid.grid_size = [nx,ny,1];
    end
    
    M_rigid = zeros([nx, ny, perframes(i)],'single');
    M_rigid = apply_shifts(M,shifts_rigid,options_rigid);
    M_rigid2 = zeros([nx2, ny2, perframes(i)],'single');
    M_rigid2 = apply_shifts(M2,shifts_rigid,options_rigid);
    clear M M2
    clear shifts_rigid

    % NON-rigid motion correction starts below -------------------------
    if nonrigid_mc
        disp(sprintf('\t \t %s: Starting non-rigid motion correction ', datestr(now)))

        M_proc = M_rigid;
        M_proc2 = M_rigid2;
        if bandpass
            M_proc = spatial_bandpass(M_proc,avg_cell_radius,10,2,use_gpu);
            M_proc2 = spatial_bandpass(M_proc2,avg_cell_radius,10,2,use_gpu);
        end

        options_nonrigid = NoRMCorreSetParms('d1',nx,'d2',ny,'max_shift',max_shift,'us_fac',30,'grid_size',[ns_nonrigid,ns_nonrigid],'print_msg',0);
        [~,shifts_nonrigid,~,options_nonrigid] = normcorre_batch(M_proc,options_nonrigid,im_ds_nonrig); 
        clear M_proc M_proc2

        M_final = zeros([nx, ny, perframes(i)],'single');
        M_final = apply_shifts(M_rigid,shifts_nonrigid,options_nonrigid);
        M_final2 = zeros([nx2, ny2, perframes(i)],'single');
        M_final2 = apply_shifts(M_rigid2,shifts_nonrigid,options_nonrigid);
        clear shifts_nonrigid
        clear M_rigid M_rigid2
    else
        % skip non-rigid if nonrigid_mc == 0
        disp(sprintf('\t \t %s: No non-rigid motion correction', datestr(now)))
        M_final= M_rigid;
        M_final2= M_rigid2;
        clear M_rigid M_rigid2
    end

    disp(sprintf('\t \t %s: Saving the motion corrected movie ', datestr(now)))

    h5write(output_filename,output_datasetname,single(M_final),[1,1,startno(i)],[nx,ny,perframes(i)]);
    h5write(output_filename2,output_datasetname2,single(M_final2),[1,1,startno(i)],[nx2,ny2,perframes(i)]);

    clear M_final M_final2
end
disp(sprintf('%s: Motion correction finished', datestr(now)))

end
   
