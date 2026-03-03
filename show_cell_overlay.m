function show_cell_overlay(output, idx, varargin)
    
    contour_thresh = 0.3;
    clim_scale = [0, 1];
    draw_hole = 0;
    label = 0;
    base = output.info.summary_image;
    ftsz = 12;
    if nargin>2
        for k = 1:length(varargin)
            vararg = varargin{k};
            if ischar(vararg)
                switch lower(vararg)
                    case 'base'
                        base = double(varargin{k+1});
                    case 'clim_scale'
                        clim_scale = varargin{k+1};
                    case 'contour_thresh'
                        contour_thresh = varargin{k+1};
                    case 'draw_hole'
                        draw_hole = varargin{k+1};
                    case 'label'
                        label = varargin{k+1};
                    case 'font_size'
                        ftsz = varargin{k+1};
                end
            end
        end
    end

    clims = quantile(base(:), clim_scale);
    imagesc(base, clims); axis image; axis off;
    colormap('gray');
    for i = 1:length(idx)
        c = max(rand(1, 3), 0.2);
        c = c / max(c);
        cell_image = full(full(output.spatial_weights(:,:,idx(i))));
        max_val = max(max(cell_image));
        b = bwboundaries(cell_image > contour_thresh * max_val,'holes');
        lens = cellfun(@length, b);
        if ~isempty(lens)
            if (~draw_hole) || length(b)<2
                b = b{find(lens == max(lens), 1)};
                center = mean(b,1);
                hold on;
                plot((b(:,2)), (b(:,1)),'Color', c, 'LineWidth', 1);
                if label
                    text(center(2),center(1),num2str(idx(i)),"FontSize",ftsz);
                end
                hold off
            else
                [~, b_sort] = sort(lens, 'descend');
                for j=1:2
                    hold on;
                    center = mean(b{b_sort(j)},1);
                    plot((b{b_sort(j)}(:,2)), (b{b_sort(j)}(:,1)),'Color', c, 'LineWidth', 1);
                    if label
                        text(center(2),center(1),num2str(idx(i)),"FontSize",ftsz);
                    end
                    hold off
                end
            end
        end


    end

end