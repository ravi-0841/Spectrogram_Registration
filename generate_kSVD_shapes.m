param = struct;
param.K = 32;
param.numIteration = 250;
param.displayProgress = 1;
param.preserveDCAtom = 0;
param.InitializationMethod = 'DataElements';
param.L = 3;
param.errorFlag = 0;

[dict_sad, ~] = KSVD(shapes_sad,param);
[dict_happy, ~] = KSVD(shapes_happy,param);
[dict_nutrl, ~] = KSVD(shapes_nutrl,param);


%% Testing coeff sanity by generating inverse Fourier co-ordinates
num_coeffs = 128;
T = 99;
s = 1;
figure();

shapes_xy_sad = zeros(2*(T+1), size(shapes_sad,2));
shapes_xy_happy = zeros(2*(T+1), size(shapes_happy,2));
shapes_xy_nutrl = zeros(2*(T+1), size(shapes_nutrl,2));

for objs = 1:size(shapes_sad,2)
    z1 = shapes_sad(:,objs);
    pts = getPointsEllipticFourierDescriptor(z1(1:num_coeffs),...
            z1(1+num_coeffs:2*num_coeffs), z1(2*num_coeffs+1: 3*num_coeffs),...
            z1(3*num_coeffs+1:4*num_coeffs), T, s);
    shapes_xy_sad(:,objs) = [real(pts) imag(pts)]';
end

for objs = 1:size(shapes_happy,2)
    z2 = shapes_happy(:,objs);
    pts = getPointsEllipticFourierDescriptor(z2(1:num_coeffs),...
            z2(1+num_coeffs:2*num_coeffs), z2(2*num_coeffs+1: 3*num_coeffs),...
            z2(3*num_coeffs+1:4*num_coeffs), T, s);
    shapes_xy_happy(:,objs) = [real(pts) imag(pts)]';
end

for objs = 1:size(shapes_nutrl,2)
    z3 = shapes_nutrl(:,objs);
    pts = getPointsEllipticFourierDescriptor(z3(1:num_coeffs),...
            z3(1+num_coeffs:2*num_coeffs), z3(2*num_coeffs+1: 3*num_coeffs),...
            z3(3*num_coeffs+1:4*num_coeffs), T, s);
    shapes_xy_nutrl(:,objs) = [real(pts) imag(pts)]';
end

%%
for i = 1:100
    z1 = shapes_sad(:,i);
    z2 = shapes_happy(:,i);
    z3 = shapes_nutrl(:,i);
    x1 = ellipticFourierFunction(z1(1:num_coeffs),...
            z1(1+num_coeffs:2*num_coeffs),100);
    y1 = ellipticFourierFunction(z1(2*num_coeffs+1: 3*num_coeffs),...
            z1(3*num_coeffs+1:4*num_coeffs),100);
    k = 1;
    pl = [];
    for it = 0:s:T
        pl(k) = x1(it) + 1i*y1(it);
        k = k + 1;
    end
    subplot(131), plot(real(pl), imag(pl)), title('Sad');
    
    x2 = ellipticFourierFunction(z2(1:num_coeffs),...
            z2(1+num_coeffs:2*num_coeffs),100);
    y2 = ellipticFourierFunction(z2(2*num_coeffs+1: 3*num_coeffs),...
            z2(3*num_coeffs+1:4*num_coeffs),100);
    k = 1;
    pl = [];
    for it = 0:s:T
        pl(k) = x2(it) + 1i*y2(it);
        k = k + 1;
    end
    subplot(132), plot(real(pl), imag(pl)), title('Happy');
    
    x3 = ellipticFourierFunction(z3(1:num_coeffs),...
            z3(1+num_coeffs:2*num_coeffs),100);
    y3 = ellipticFourierFunction(z3(2*num_coeffs+1: 3*num_coeffs),...
            z3(3*num_coeffs+1:4*num_coeffs),100);
    k = 1;
    pl = [];
    for it = 0:s:T
        pl(k) = x3(it) + 1i*y3(it);
        k = k + 1;
    end
    subplot(133), plot(real(pl), imag(pl)), title('Neutral');
    pause;
end

%% plot x-y shapes
figure();
for i = 1:100
    subplot(131), plot(shapes_xy_sad(1:(T+1),i), shapes_xy_sad((T+2):end,i)), ...
        title('sad');
    subplot(132), plot(shapes_xy_happy(1:(T+1),i), shapes_xy_happy((T+2):end,i)), ...
        title('happy');
    subplot(133), plot(shapes_xy_nutrl(1:(T+1),i), shapes_xy_nutrl((T+2):end,i)), ...
        title('neutral');
    pause;
end