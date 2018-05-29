%% Register two images
function [Mp,sx,sy,vx,vy] = register(F,M,opt)

    if nargin<3;  opt = struct();  end;
    if ~isfield(opt,'sigma_fluid');      opt.sigma_fluid     = 1.0;              end;
    if ~isfield(opt,'sigma_diffusion');  opt.sigma_diffusion = 1.0;              end;
    if ~isfield(opt,'sigma_i');          opt.sigma_i         = 1.0;              end;
    if ~isfield(opt,'sigma_x');          opt.sigma_x         = 1.0;              end;
    if ~isfield(opt,'niter');            opt.niter           = 250;              end;
    if ~isfield(opt,'vx');               opt.vx              = zeros(size(M));   end;
    if ~isfield(opt,'vy');               opt.vy              = zeros(size(M));   end;
    if ~isfield(opt,'stop_criterium');   opt.stop_criterium  = 0.01;             end;
    if ~isfield(opt,'imagepad');         opt.imagepad        = 1.2;              end;
    if ~isfield(opt,'do_display');       opt.do_display      = 1;                end;
    if ~isfield(opt,'do_plotenergy');    opt.do_plotenergy   = 1;                end;

    %% padded image
    [F,lim] = imagepad(F,opt.imagepad);
    [M,lim] = imagepad(M,opt.imagepad);
    
    %% T is the deformation from M to F
    vx = imagepad(opt.vx,opt.imagepad);
    vy = imagepad(opt.vy,opt.imagepad);
    e  = zeros(1,opt.niter);
    e_min = 1e+100;      % Minimal energy
    
    %% Iterate update fields
    for iter=1:opt.niter

        % Find update
        [ux,uy] = findupdate(F,M,vx,vy,opt.sigma_i,opt.sigma_x);

        % Regularize update
        ux    = imgaussian(ux,opt.sigma_fluid);
        uy    = imgaussian(uy,opt.sigma_fluid);

        % Compute step (e.g., max half a pixel)
        step  = opt.sigma_x;
        
        % Update velocities (demons) - additive
        vx = vx + step*ux;
        vy = vy + step*uy;

        % Update velocities (demons) - composition
        %[vx,vy] = compose(vx,vy,step*ux,step*uy);
        
        % Regularize velocities
        vx = imgaussian(vx,opt.sigma_diffusion);
        vy = imgaussian(vy,opt.sigma_diffusion);
        
        % Get Transformation
        [sx,sy] = expfield(vx,vy);  % deformation field

        % Compute energy
        e(iter) = energy(F,M,sx,sy,opt.sigma_i,opt.sigma_x);
        disp(['Iteration: ' num2str(iter) ' - ' 'energy: ' num2str(e(iter))]);
        if e(iter)<e_min
            sx_min = sx; sy_min = sy; % update best fields
            vx_min = vx; vy_min = vy; % update best fields
            e_min  = e(iter);
        end
        
        % Stop criterium
        if iter>1 && abs(e(iter) - e(max(1,iter-5))) < e(1)*opt.stop_criterium
            break;
        end

        if opt.do_display
            % display deformation
            subplot(2,4,7); showvector(ux,uy,4,3,lim); title('Update');
            subplot(2,4,8); showgrid  (sx,sy,4,lim); title('Transformation');
            drawnow;
            
            % Display registration
            Mp     = iminterpolate(M,sx,sy);
            diff   = (F-Mp).^2;
            showimage(F,'Fixed', M,'Moving', Mp,'Warped', diff,'Diff', 'lim',lim,'nbrows',2); drawnow;

            % Plot energy
            if opt.do_plotenergy
                subplot(2,2,3)
                hold on;
                plot(1:iter,e(1:iter),'r-'); xlim([0 opt.niter]);
                xlabel('Iteration'); ylabel('Energy');
                hold off;
                drawnow
            end
        end

    end
    
    %% Get Best Transformation
    vx = vx_min;  vy = vy_min;
    sx = sx_min;  sy = sy_min;
    
    %% Transform moving image
    Mp = iminterpolate(M,sx,sy);
    
    %% Unpad image
    Mp = Mp(lim(1):lim(2),lim(3):lim(4));
    vx = vx(lim(1):lim(2),lim(3):lim(4));
    vy = vy(lim(1):lim(2),lim(3):lim(4));
    sx = sx(lim(1):lim(2),lim(3):lim(4));
    sy = sy(lim(1):lim(2),lim(3):lim(4));

end