function [traces] = getRawTraces(movie_in, filters_in, varargin)
% Extract raw signals by straightforward multiplication
% between raw image data and spatial filter

% ZW 20260212

    if ischar(movie_in)
        movie_dims = h5info(movie_in).Datasets.Dataspace.Size;
        get_frames = @(x) h5read(movie_in, '/mov', [1 1 x(1)], [movie_dims(1) movie_dims(2) x(end)-x(1)+1]);
    else
        movie_dims = size(movie_in);
        get_frames = @(x) movie_in(:,:,x);
    end
    height = movie_dims(1);
    width = movie_dims(2);
    num_pixels = height * width;
    num_frames = movie_dims(3);

    % figure; imagesc(get_frames(1))

    get_filter = @(x) filters_in(:,:,x);  
    num_filters = size(filters_in,3);
    filters = zeros(num_pixels, num_filters, 'single'); % Preallocate
    for k = 1:num_filters
        filter = get_filter(k);
        filters(:,k) = reshape(filter, num_pixels, 1);
    end
    
    traces = zeros(num_filters, num_frames, 'single');

    frame_chunk_size = 1000;
    [frame_chunks, num_chunks] = make_frame_chunks(num_frames, frame_chunk_size);


    for i = 1:num_chunks
        frame_inds = frame_chunks(i,1):frame_chunks(i,2); % For this chunk
        num_frames_in_chunk = frame_inds(end) - frame_inds(1) + 1;
        fprintf('%s: Calculating traces for frames %d to %d (out of %d) \n',...
            datestr(now), frame_inds(1), frame_inds(end), num_frames);
            
        movie_chunk = get_frames(frame_inds);
        movie_chunk = reshape(movie_chunk, num_pixels, num_frames_in_chunk);

        traces(:,frame_inds) = filters' * single(movie_chunk);
    end

    traces = traces'; % [num_frames x num_filters]

end