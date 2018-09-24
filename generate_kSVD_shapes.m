param = struct;
param.K = 32;
param.numIteration = 200;
param.displayProgress = 1;
param.preserveDCAtom = 0;
param.InitializationMethod = 'DataElements';
param.L = 3;
param.errorFlag = 0;

[dict_sad, ~] = KSVD(shapes_sad,param);
[dict_happy, ~] = KSVD(shapes_happy,param);
[dict_nutrl, ~] = KSVD(shapes_nutrl,param);