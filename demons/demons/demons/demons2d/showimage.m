%% Display two images
%  Changed: Dec 6th, 2011
%
function showimage(varargin)

    % Check parameters
    nb_args   = size(varargin,2);
    nb_images = nb_args;
    nb_rows   = 1;
    row       = 1;
    crange    = [0 1]*256; % default image intensities
    
    for i=1:nb_args
        if ischar(varargin{i})
            if isequal(varargin{i},'lim')
                lim       = varargin{i+1};
                nb_images = nb_images-2;
            elseif isequal(varargin{i},'nbrows')
                nb_rows   = varargin{i+1};
                nb_images = nb_images-2;
            elseif isequal(varargin{i},'row')
                row       = varargin{i+1};
                if row>nb_rows; nb_rows = row; end;
                nb_images = nb_images-2;
            elseif isequal(varargin{i},'caxis')
                crange    = varargin{i+1};
                nb_images = nb_images-2;
            else
                nb_images = nb_images-1;
            end
        end
    end
    
    % Display images
    iter_image = 1;
    for iter_arg=1:nb_args
        if ~ischar(varargin{iter_arg})
            I = varargin{iter_arg};
            subplot(nb_rows,nb_images,(row-1)*nb_images + iter_image);
            imagesc(I,crange);
            daspect([1 1 1]);
            if exist('lim'); axis([lim(3) lim(4) lim(1) lim(2)]); end
            axis off;
            if iter_arg+1<=nb_args && ischar(varargin{iter_arg+1})
                title(varargin{iter_arg+1});
            end
            iter_image = iter_image+1;
        end
        if iter_image>nb_images
            break;
        end
    end
        
end