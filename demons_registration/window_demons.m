function disp_field = window_demons(F, M, opts)
    if nargin<3;                         opts               = struct();     end
    if ~isfield(opts,'window_size');     opts.window        = 4;            end
    if ~isfield(opts,'stride');          opts.stride        = 2;            end
    if ~isfield(opts,'alpha');           opts.alpha         = 0.4;          end
    if ~isfield(opts,'sigma_fluid');     opts.sigma_fluid   = 1.0;          end
    if ~isfield(opts,'sigma_diff');      opts.sigma_diff    = 1.0;          end
    if ~isfield(opts,'lambda');          opts.lambda        = 0.1;          end
    if ~isfield(opts,'step');            opts.step          = 1.0;          end
    if ~isfield(opts,'epsilon');         opts.epsilon       = 10;           end
    if ~isfield(opts,'max_iter');        opts.max_iter      = 10;           end
    if ~isfield(opts,'max_epochs');      opts.max_epochs    = 3;            end

    org_size = size(F,2);

    num_cols    = size(F, 2);
    num_frames  = ceil((num_cols - opts.window)/opts.stride) + 1;
    extend_cols = (num_frames - 1)*opts.stride + opts.window;
    padding     = extend_cols - num_cols;

    F = [F zeros(size(F,1),padding)];
    M = [M zeros(size(M,1),padding)];

    fdf = cell(num_frames,2);

    for frame = 1:num_frames
        fdf{frame,1} = zeros(size(F));
        fdf{frame,2} = zeros(size(F));
    end

    sdf = cell(num_frames,2);

    for frame = 1:num_frames
        sdf{frame,1} = zeros(size(F));
        sdf{frame,2} = zeros(size(F));
    end

    epochs = 1;
    while (epochs < opts.max_epochs)
        cur_frame = 1;
        disp(['Current Iteration: ' num2str(epochs)]);

        for start_index = 1:opts.stride:(extend_cols - opts.stride)
            
            F_tilda = zeros(size(F));
            F_tilda(:, start_index:start_index+opts.window-1) = ...
                        F(:, start_index:start_index+opts.window-1);
            
            M_tilda = zeros(size(M));
            M_tilda(:, start_index:start_index+opts.window-1) = ...
                        M(:, start_index:start_index+opts.window-1);
            
            [G_fix_x, G_fix_y] = imgradientxy(F_tilda, 'central');
            [G_fix_mag, ~] = imgradient(G_fix_x, G_fix_y);

            num_iter = 1;

            current_moved = imwarp(M_tilda, ...
                cat(3,fdf{cur_frame,1},fdf{cur_frame,2}));

            while (num_iter < opts.max_iter)

                [G_mov_x, G_mov_y] = imgradientxy(current_moved, 'central');
                [G_mov_mag, ~] = imgradient(G_mov_x, G_mov_y);

                I_diff = current_moved - F_tilda;
                normzf = opts.alpha.^2 * I_diff.^2 + G_fix_mag.^2 + opts.lambda;
                normzm = opts.alpha.^2 * I_diff.^2 + G_mov_mag.^2 + opts.lambda;

                if cur_frame>1
                    u_x = -1 * ((I_diff.*G_fix_x ...
                        + opts.lambda*sdf{cur_frame-1,1})./normzf ...
                        + (I_diff.*G_mov_x ...
                        + opts.lambda*sdf{cur_frame-1,1})./normzm);

                    u_y = -1 * ((I_diff.*G_fix_y ...
                        + opts.lambda*sdf{cur_frame-1,2})./normzf ...
                        + (I_diff.*G_mov_y ...
                        + opts.lambda*sdf{cur_frame-1,2})./normzm); 

                else
                    u_x = -1 * ((I_diff.*G_fix_x)./normzf ...
                        + (I_diff.*G_mov_x)./normzm);

                    u_y = -1 * ((I_diff.*G_fix_y)./normzm ...
                        + (I_diff.*G_mov_y)./normzm);
                end

                u_x(isnan(u_x)) = 0;
                u_y(isnan(u_y)) = 0;

                sdf{cur_frame,1} = imgaussfilt(u_x, opts.sigma_fluid);
                sdf{cur_frame,2} = imgaussfilt(u_y, opts.sigma_fluid);

                fdf{cur_frame,1} = imgaussfilt(fdf{cur_frame,1} + u_x, ...
                                    opts.sigma_diff);
                fdf{cur_frame,2} = imgaussfilt(fdf{cur_frame,2} + u_y, ...
                                    opts.sigma_diff);

                current_moved = imwarp(M_tilda, ...
                                cat(3,fdf{cur_frame,1},fdf{cur_frame,2}));

                num_iter = num_iter + 1;
            end
            cur_frame = cur_frame + 1;
        end
        epochs = epochs + 1;
    end

    stitch_dfx = zeros(size(F));
    stitch_dfy = zeros(size(F));
    counter = 1;
    for start_index = 1:opts.stride:(extend_cols - opts.stride)
        stitch_dfx(:,start_index:start_index+opts.window-1) = ...
        fdf{counter,1}(:,start_index:start_index+opts.window-1);
        stitch_dfy(:,start_index:start_index+opts.window-1) = ...
        fdf{counter,2}(:,start_index:start_index+opts.window-1);
        counter = counter + 1;
    end
    stitch_dfx = imgaussfilt(stitch_dfx,opts.sigma_diff);
    stitch_dfy = imgaussfilt(stitch_dfy,opts.sigma_diff);
    disp_field = cat(3,stitch_dfx(:,1:org_size),stitch_dfy(:,1:org_size));
    figure(), showvector(squeeze(disp_field(:,:,1)),squeeze(disp_field(:,:,2)),5);
end