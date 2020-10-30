%% Home Assignment 2 - FFR135 Artificial Neural Networks
% Author: Arshad Nowsath(nowsath)
%% Loading the dataset from MNIST

clc; 
clear;
clear all;

[xTrain, tTrain, xValid, tValid, xTest, tTest] = LoadMNIST(3);

%% Convolutional neural network 1

% Options to train the network using stochastic gradient descent

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

% Layout of the layers
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

% Training the network
network_1 = trainNetwork(xTrain, tTrain, layers, options);

% Computing the scores
[prediction_train_1, scores_train_1] = classify(network_1, xTrain);
[prediction_valid_1, scores_valid_1] = classify(network_1, xValid);
[prediction_test_1, scores_test_1] = classify(network_1, xTest);

% Classification errors  obtained on the training, validation, and test sets sets
Classification_train_1 = classification_error(tTrain, prediction_train_1);
Classification_valid_1 = classification_error(tValid, prediction_valid_1);
Classification_test_1 = classification_error(tTest, prediction_test_1);

%% Convolutional neural network 1

% Options to train the network using stochastic gradient descent

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

% Layout of the layers
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

% Training the network
network_2 = trainNetwork(xTrain, tTrain, layers, options);

% Computing the scores
[prediction_train_2, scores_train_2] = classify(network_2, xTrain);
[prediction_valid_2, scores_valid_2] = classify(network_2, xValid);
[prediction_test_2, scores_test_2] = classify(network_2, xTest);

% Classification errors  obtained on the training, validation, and test sets sets
Classification_train_2 = classification_error(tTrain, prediction_train_2);
Classification_valid_2 = classification_error(tValid, prediction_valid_2);
Classification_test_2 = classification_error(tTest, prediction_test_2);

%% Saving the network

save network_1;

save network_2;
%% Classification error function

function C = classification_error(target, outputs)
length_valset = size(target,1);    
C = 1 / length_valset * sum(outputs ~= target);
end
