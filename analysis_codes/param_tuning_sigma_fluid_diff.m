%% SSD based registration Parameter Tuning

warning('off', 'all');
[s_ang,~] = audioread('./wav/angry.wav');
[s_neu,~] = audioread('./wav/neutral.wav');

f = 16000;
r = 512;
w = 0.025;
s = 0.015;
linear_registration = 0;

[s_neu, s_ang] = get_alignment(s_neu,s_ang,f,w,w-s,r,3);

spect_ang = spectrogram(s_ang,w*f,int64((w-s)*f),r);
spect_neu = spectrogram(s_neu,w*f,int64((w-s)*f),r);

spect_ang_mag = abs(spect_ang);
spect_ang_phase = angle(spect_ang);

spect_neu_mag = abs(spect_neu);
spect_neu_phase = angle(spect_neu);

spect_neu_mag = mat2gray(spect_neu_mag);
spect_ang_mag = mat2gray(spect_ang_mag);


sigma_diff_grid = 0.5:0.1:2;
sigma_fluid_grid = 0.5:0.1:2;

ssd_matrix = zeros(length(sigma_diff_grid), length(sigma_fluid_grid));
mi_matrix = zeros(length(sigma_diff_grid), length(sigma_fluid_grid));

for sd = 1:length(sigma_diff_grid)
    for sf = 1:length(sigma_fluid_grid)
        disp(['Sigma Fluid: ' num2str(sigma_fluid_grid(sf)) ' and ' 'Sigma Diff: ' num2str(sigma_diff_grid(sd))]);
        [~,~,ssd_matrix(sd,sf),mi_matrix(sd,sf)] = my_demons(spect_ang_mag, spect_neu_mag, ...
                                    0.4, sigma_fluid_grid(sf), sigma_diff_grid(sd), 0.000001);
    end
end