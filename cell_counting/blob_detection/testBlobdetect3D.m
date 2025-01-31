
function tests = testBlobdetect3D
tests = functiontests(localfunctions);
end

%% argparser
function testArgs(testCase)
%takes reasonable values for all keyword arguments in doc
array = zeros([20,20,20],'uint8');
blobs=blobdetect3D(array, 12, Filter='LoG', DarkBackground=true, MedianFilter=false, ...
    QualityFilter=0.3, OverlapFilter = true, KernelSize=6, GPU=true);
end

function test2DArray(testCase)
%error in case of 2D input
array = zeros([20,20],'uint8');
verifyError(testCase, @() blobdetect3D(array, 12), "MATLAB:InputParser:ArgumentFailedValidation");
end

%% GPU and GPU memory
function testSmallData(testCase)
smalldata = zeros([20,20,20],'uint8');
verifyWarningFree(testCase, @() blobdetect3D(smalldata,16, GPU=true));
end

function testSmallDataNoGPU(testCase)
smalldata = zeros([20,20,20],'uint8');
%default no GPU, convn used
verifyWarningFree(testCase, @() blobdetect3D(smalldata,16));
end

function testBigData(testCase)
bigdata = zeros([5000,2000,200],'uint8');
verifyError(testCase, @() blobdetect3D(bigdata,16, GPU=true), "parallel:gpu:array:OOM");
end

% function testBigDataNoPGU(testCase)
% bigdata = zeros([5000,2000,200],'uint8');
% takes too long to test, ~40min
% verifyWarningFree(testCase, @() blobdetect3D(bigdata,16));
% end

%% output format

function testOutputFormat(testCase)
array = zeros([20,20,20],'uint8');
blobs=blobdetect3D(array, 12);
%check that a table (probably empty) with headers x y z Intensity came out
assert(all(ismember({'x','y','z','Intensity'}, blobs.Properties.VariableNames)))
end

%% find all blobs in simple test data
