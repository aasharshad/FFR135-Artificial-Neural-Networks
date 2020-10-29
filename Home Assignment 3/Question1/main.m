clc; clear all;
[xTrain, tTrain, xValid, tValid, xTest, tTest] = LoadMNIST(3);

%tTrain=tTrain.';

options = trainingOptions('sgdm', ...
    'Momentum', 0.9, ...
    'MaxEpochs',60, ...    
    'InitialLearnRate',0.001, ...
    'MiniBatchSize',8192, ...
    'Shuffle', 'every-epoch', ...
    'ValidationData',{xValid,tValid}, ...
    'ValidationFrequency',30, ...
    'ValidationPatience',5, ...
    'Plots','training-progress');

% Define layers
layers = [
    imageInputLayer([28 28 1])
    
    convolution2dLayer(5,20,'stride',1,'Padding',1,'WeightsInitializer','narrow-normal')
    
    reluLayer
    
    maxPooling2dLayer(2,'Stride',2) 
    
    fullyConnectedLayer(100,'WeightsInitializer','narrow-normal')
    
    reluLayer
    
    fullyConnectedLayer(10,'WeightsInitializer','narrow-normal')
    softmaxLayer
    
    classificationLayer];

% Train network
net_1 = trainNetwork(xTrain, tTrain, layers, options);

% Compute scores
[pred_train_1, scores_train_1] = classify(net_1, xTrain);
[pred_valid_1, scores_valid_1] = classify(net_1, xValid);
[pred_test_1, scores_test_1] = classify(net_1, xTest);

C_train_1 = classification_error(tTrain, pred_train_1);
C_valid_1 = classification_error(tValid, pred_valid_1);
C_test_1 = classification_error(tTest, pred_test_1);

%%

options = trainingOptions('sgdm', ...
    'Momentum', 0.9, ...
    'InitialLearnRate',0.01, ...
    'MaxEpochs',30, ...
    'MiniBatchSize',8192, ...
    'Shuffle', 'every-epoch', ...
    'ValidationData',{xValid,tValid}, ...
    'ValidationFrequency',30, ...
    'ValidationPatience',5, ...
    'Plots','training-progress');

% Define layers
layers = [
    imageInputLayer([28 28 1])
    
    convolution2dLayer(3,20,'stride',1,'Padding',1,'WeightsInitializer','narrow-normal')
    batchNormalizationLayer
    reluLayer
    
    maxPooling2dLayer(2,'Stride',2) 
    
    convolution2dLayer(3,30,'stride',1,'Padding',1,'WeightsInitializer','narrow-normal')
    batchNormalizationLayer
    reluLayer
    
    maxPooling2dLayer(2,'Stride',2) 
    
    convolution2dLayer(3,50,'stride',1,'Padding',1,'WeightsInitializer','narrow-normal')
    batchNormalizationLayer
    reluLayer
      
    
    fullyConnectedLayer(10,'WeightsInitializer','narrow-normal')
    softmaxLayer
    
    classificationLayer];

% Train network
net_2 = trainNetwork(xTrain, tTrain, layers, options);

% Compute scores
[pred_train_2, scores_train_2] = classify(net_2, xTrain);
[pred_valid_2, scores_valid_2] = classify(net_2, xValid);
[pred_test_2, scores_test_2] = classify(net_2, xTest);

C_train_2 = classification_error(tTrain, pred_train_2);
C_valid_2 = classification_error(tValid, pred_valid_2);
C_test_2 = classification_error(tTest, pred_test_2);

%%
save net_2;

save net_1;