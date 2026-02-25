function save_h5_range(h5_in, range, h5_out)


    chunk_size=1000;
    
    
    %%
    movie_dims = h5info(h5_in).Datasets.Dataspace.Size;
    height = movie_dims(1);
    width = movie_dims(2);
    chunk_first_last = chunkFrames(chunk_size,[range(1) range(2)]);
    
    get_frames = @(x) h5read(h5_in, '/mov', [1 1 x(1)], [height width x(end)-x(1)+1]);
    if isfile(h5_out)
        disps(['Already found h5 file:' h5_out 'deleting!']);
        delete(h5_out);
    end
    disp('Saving part of h5 file...')
    
    for ichunk=1:size(chunk_first_last,1) 
        movieBatch = get_frames([chunk_first_last(ichunk,1),chunk_first_last(ichunk,2)]);
        h5append(h5_out, single(movieBatch)); 
    end
    
    disp('Saving done.')

end