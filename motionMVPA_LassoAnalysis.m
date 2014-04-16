motionMVPA_CreateMetadata
addpath('~/matlab/glmnet');

%% Load Data
X = zeros(2880,55698);
for ss=1:12
	a=((ss-1)*240)+1;
	b=240*ss;
	load(sprintf('tSeries%d.mat',ss));
	X(a:b,:) = bsxfun(@minus,tSeries,mean(tSeries));
end
fprintf('Data loaded\n');

%% Generate options structures
optsCV = glmnetSet();
optsFit = glmnetSet();

%% HRF Shift/Adjustments
% temp = reshape(y,240,12);
% temp(4:end,:) = temp(1:end-3,:);
% temp(1:3,:) = 0;
% y = temp(:);

% Average 3rd and 4th TR after Stimulus Onset
%        v v
% SO 1 2 3 4 (TR/conceptual idexing)
%  1 2 3 4 5 (Matlab indexing)
%
% 240/3 = 80 stimulus presentations
% It is certain that the last stimulus presentation will not have data
% under this method of aggregation.
X = mat2cell(X,repmat(240,12,1),size(X,2));

ix = bsxfun(@plus,[4;5],3*(0:(80-2))); % -1 for zero base, -1 to drop last trial.


X = cellfun(@(x) x(ix(:),:),X,'unif',0);
X = cell2mat(X);
X = cellfun(@mean,mat2cell(X,repmat(2,(80-1)*12,1),size(X,2)),'unif',0);
X = cell2mat(X);

z = true(2880,1);
z = z & (mod(1:2880,3)==1)';
z = z & (mod(0:2880-1,240)<237)';
z = z & metadata.True2D;

%% Apply data reduction/filtering
% Only 2d motion
X = X(metadata.True2D(z),:);
y = double(metadata.TrueLeft(z));

%% Make CV
[~,cv] = ndgrid(1:240,[1:12]);
cv = cv(:);
cv = cv(z);
train = cv < 12;
test = cv == 12;
cv = cv(cv<12);

%% Run 
Xtrain = X(train,1:1000);
ytrain = y(train);
glmnetCV = cvglmnet(Xtrain,ytrain,'binomial',optsCV,'class',11,cv',false,false);
optsFit.lambda = glmnetCV.lambda_min;
glmnetFit = glmnet(Xtrain,ytrain,'binomial',optsFit);
