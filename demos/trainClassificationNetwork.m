
% the text file contains paths of training data produced by CellLabelingTool
textfilewithpaths = 'D:\DATA_folder\training_datasets.txt';
netowrksavepath   = 'D:\DATA_folder\'; % where to save the trained network
%==========================================================================
dp = importdata(textfilewithpaths);

XTrain = cell(numel(dp), 1);
YTrain = cell(numel(dp), 1);
for ii = 1:numel(dp)
    dataload = load(dp{ii});
    XTrain{ii} = dataload.XTrain;
    YTrain{ii} = dataload.labels;
end
XTrainall = cat(4, XTrain{:});
YTrainall = cat(1, YTrain{:});
idxuse = balanceChoices(YTrainall);
XTrain = XTrainall(:, :, :, idxuse);
YTrain = categorical(YTrainall(idxuse));

imgsize = size(XTrain, 1:3);
% --- 1. Split Data into Train and Validation (80/20 split) ---
idx = randperm(size(XTrain, 4));
nTrain = floor(0.8 * length(idx));

idxTrain = idx(1:nTrain);
idxVal   = idx(nTrain+1:end);

X_tr = XTrain(:, :, :, idxTrain);
Y_tr = YTrain(idxTrain);
X_val = XTrain(:, :, :, idxVal);
Y_val = YTrain(idxVal);

% --- 2. Define Augmentation ---
% Adds rotation, reflection, and translation to make the net robust
dsX = arrayDatastore(X_tr, 'IterationDimension', 4);
dsY = arrayDatastore(Y_tr);
dsTrainRaw = combine(dsX, dsY);
augimdsTrain = transform(dsTrainRaw, @augmentCellViews);
% augmenter = imageDataAugmenter( ...
%     'RandRotation', [-20 20], ...      % Rotate +/- 20 degrees
%     'RandScale', [0.8 1.2], ...
%     'RandXReflection', true, ...       % Flip horizontally
%     'RandYReflection', true); ...       % Flip vertically);
% 
% % Create a datastore that applies these augmentations on the fly
% augimdsTrain = augmentedImageDatastore(imgsize, X_tr, Y_tr, ...
%     'DataAugmentation', augmenter);

% --- 3. Define Architecture ---
layers = [
    imageInputLayer(imgsize) % 3 channels for the 3 views
    
    convolution2dLayer(3, 16, 'Padding', 'same')
    batchNormalizationLayer
    reluLayer
    maxPooling2dLayer(2, 'Stride', 2)
    
    convolution2dLayer(3, 32, 'Padding', 'same')
    batchNormalizationLayer
    reluLayer
    maxPooling2dLayer(2, 'Stride', 2)
    
    convolution2dLayer(3, 64, 'Padding', 'same')
    batchNormalizationLayer
    reluLayer
    
    fullyConnectedLayer(2)
    softmaxLayer
    classificationLayer
];

% --- 4. Train ---
options = trainingOptions('adam', ...
    'MaxEpochs', 100, ...
    'MiniBatchSize', 64, ...
    'ValidationData', {X_val, Y_val}, ... % Validate on non-augmented data
    'Plots', 'training-progress', ...
    'Shuffle', 'every-epoch', ...
    'Verbose', false);
fprintf('Training with %d samples...\n', numel(idx))
%%
net    = trainNetwork(augimdsTrain, layers, options);
% net    = trainNetwork(X_tr, Y_tr, layers, options);

% --- 5. Test ---
YPred       = classify(net, X_val);
valaccuracy = mean(YPred == Y_val);
YPredtr     = classify(net, X_tr);
traccuracy  = mean(YPredtr == Y_tr);

YPredall    = classify(net, XTrainall);
totaccuracy = mean(YPredall == categorical(YTrainall));


save(fullfile(netowrksavepath,...
    sprintf('%s_CellClassifierNet.mat', string(datetime('now'),'yyyyMMdd'))), ...
    'net', 'valaccuracy', 'totaccuracy', 'traccuracy');
fprintf('Training succeeded! Train error: %2.3f. Test error: %2.3f\n', traccuracy, valaccuracy)

visualizeNetworkClassification(net, X_val, Y_val, netowrksavepath);

