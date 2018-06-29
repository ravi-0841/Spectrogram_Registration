If = zeros(200,200);
Im = zeros(200,200);

If(80:120,80:120) = 1;
Im(60:110,70:105) = 1;

sigma_diff_array = 1.0:0.2:2.5;
sigma_flud_array = 1.0:0.2:2.5;

ssd_dem = zeros(length(sigma_diff_array), length(sigma_diff_array));
ssd_win_dem = zeros(length(sigma_diff_array), length(sigma_diff_array));

opts = struct();
opts.alpha = 0.4;
opts.only_freq = 0;
opts.sigma_fluid = 1.0; %1.5
opts.sigma_diff = 1.0;  %2.5
opts.window_size = 100;
opts.stride = 100;
opts.max_epochs = 1;
opts.lambda = 0.07;
opts.step = 1.0;
opts.max_iter = 500;
opts.pyramid_levels  = 1;
opts.compositive = 0;
opts.diffeomorphism = 0;
opts.plot = 1;

for sd = 1:length(sigma_diff_array)
    for sf = 1:length(sigma_flud_array)
        disp(['Current sigma_diff is ', num2str(sigma_diff_array(sd)) ...
            ' and sigma_fluid is ', num2str(sigma_flud_array(sf))]);
        opts.sigma_diff = sigma_diff_array(sd);
        opts.sigma_fluid = sigma_flud_array(sf);
        opts.lambda = 0;
        disp_dem = my_multires_demons(If,Im,200,opts);
        Imm = imwarp(Im, disp_dem);
        ssd_dem(sd,sf) = sum(sum((If - Imm).^2));
        opts.lambda = 0.01;
        disp_win = my_multires_demons(If,Im,200,opts);
        Imw = imwarp(Im, disp_win);
        ssd_win_dem(sd,sf) = sum(sum((If - Imw).^2));
        disp(['Normal SSD ', num2str(ssd_dem(sd,sf)), ' Window SSD ', ...
            num2str(ssd_win_dem(sd,sf))]);
        break;
    end
    break;
end     