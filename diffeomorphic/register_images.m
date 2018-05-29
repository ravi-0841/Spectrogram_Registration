function register_images(F, M, vx, vy, opts)
    stop_flag = 0;
    n_iter = 1;
    
    while stop_flag~=1
       [ux,uy] = update_vec_field(F,M,vx,vy,opts);
       if opts.compositive
           [ux,uy] = exponential_mapping(ux,uy);
           [vx,vy] = compose_vec_fields(vx,vy,opts.step*ux,opts.step*uy);
       else
           vx = vx + opts.step*ux;
           vy = vy + opts.step*uy;
       end
       
       M_tilda = iminterpolate(M,vx,vy);
       
       figure(1)
       subplot(121), imshowpair(F, M_tilda);
       ssd = sum(sum((F - M_tilda).^2));
       subplot(122), plot(ssd);
       pause(0.001);
       
       if ssd<opts.stop_criterion || n_iter>opts.maxIter
           stop_flag = 1;
       end
       
       n_iter = n_iter + 1;
    end
end