function [vx,vy] = register_images(F, M, vx, vy, opts)
    stop_flag = 0;
    n_iter = 1;
    ssd = [];
    mi = [];
    
    figure(1)
    subplot(121), title('Image overlay');
    subplot(122), title('SSD');
    
    while stop_flag~=1
        [ux,uy] = update_vec_field(F,M,vx,vy,opts);
        
        ux = imgaussfilt(ux,opts.sigma_fluid);
        uy = imgaussfilt(uy,opts.sigma_fluid);
        
        if opts.compositive
            [ux,uy] = exponential_mapping(ux,uy);
            [vx,vy] = compose_vec_fields(vx,vy,opts.step*ux,opts.step*uy);
        else
            vx = vx + opts.step*ux;
            vy = vy + opts.step*uy;
        end
        
        vx = imgaussfilt(vx,opts.sigma_diffusion);
        vy = imgaussfilt(vy,opts.sigma_diffusion);
       
        M_tilda = iminterpolate(M,vx,vy);
%         M_tilda = imwarp(M, cat(3, vx, vy));
       
        subplot(121), imshowpair(F, M_tilda);
        ssd = [ssd sum(sum((F - M_tilda).^2))];
        mi = [mi mutual_info(F,M_tilda)];
        subplot(122), hold on;
        graph = plot(mi, 'r');
        set(graph,'LineWidth',2);
        hold off;
        pause(0.001);
       
        if ssd(end)<opts.stop_criterion || n_iter>opts.max_iter
            stop_flag = 1;
        end
       
        n_iter = n_iter + 1;
    end
end