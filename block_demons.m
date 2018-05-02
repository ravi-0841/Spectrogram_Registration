function disp_field = block_demons(moving,fixed,block,overlap,smoothing)
    assert(size(moving,1)==size(fixed,1), 'Please make sure that the moving and fixed images are same size');
    assert(size(moving,2)==size(fixed,2), 'Please make sure that the moving and fixed images are same size');
    disp_field = cat(3,zeros(size(moving)),zeros(size(moving)));
    edge_map = edge(moving,'canny');
    for i = 1:block-overlap:size(moving,1)-block
        for j = 1:block-overlap:size(moving,2)-block
            moving_block = moving(i:i+block-1,j:j+block-1);
            fixed_block = fixed(i:i+block-1,j:j+block-1);
            if sum(sum(edge_map(i:i+block-1,j:j+block-1)))>15
                [mini_disp_field,~] = imregdemons(moving_block,fixed_block,[300],'AccumulatedFieldSmoothing',2.0);
                disp_field(i:i+block-1,j:j+block-1,1) = mini_disp_field(:,:,1);
                disp_field(i:i+block-1,j:j+block-1,2) = mini_disp_field(:,:,2);
            end
        end
    end
    u = squeeze(disp_field(:,:,1));
    v = squeeze(disp_field(:,:,2));
    h = fspecial('gaussian',block,smoothing);
    u_tilda = filter2(h,u);
    v_tilda = filter2(h,v);
    disp_field(:,:,1) = u_tilda;
    disp_field(:,:,2) = v_tilda;
    abs_disp_field = squeeze(disp_field(:,:,1)).^2 + squeeze(disp_field(:,:,2)).^2;
    flow = opticalFlow(-1*squeeze(disp_field(:,:,1)),-1*squeeze(disp_field(:,:,2)));
    figure(), subplot(121), plot(flow, 'DecimationFactor',[4,16]), subplot(122), imshow(abs_disp_field), colormap('jet')
end