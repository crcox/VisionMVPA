% trial codes:
%
%   21, 22, 25, 26, 37, 38, 41, 42
%
%            L   L   R   R
%            S   F   S   F
%     2D = [22, 26, 38, 42]
%     3D = [21, 25, 37, 41]
%
%            L   L   R   R
%            3   2   3   2
%   slow = [21, 22, 37, 38]
%   fast = [25, 26, 41, 42]
%
%            2   2   3   3        
%            S   F   S   F 
%   left = [22, 26, 21, 25]
%  right = [38, 42, 37, 41]
%
%   The first 4 are dummy trials so those can go.
%
%   84 trial codes
%   Drop the first 4
%   3TRs per trial
%   240 TRs per session.
%   12 sessions * 240 TRs = 2880 TRs total;

%% Load vector data into logical TR labels 
vector_path = @(ii) sprintf('trial_vector-%d.mat',ii);
metadata.True2D = false(2880,1);
metadata.True3D = false(2880,1);
metadata.TrueSlow = false(2880,1);
metadata.TrueFast = false(2880,1);
metadata.TrueLeft = false(2880,1);
metadata.TrueRight = false(2880,1);
metadata.mriRuns = zeros(2880,1);
for ii = 1:12
    a = ((ii-1)*240)+1;
    b = ii*240;
    
    load(vector_path(ii), 'trial_vector')
    trial_vector = trial_vector(5:84);
    
    temp = any(bsxfun(@eq,trial_vector, [22, 26, 38, 42]),2);
    temp = reshape(repmat(temp',3,1),240,1);
    metadata.True2D(a:b) = temp;
    
    temp = any(bsxfun(@eq,trial_vector, [21, 25, 37, 41]),2);
    temp = reshape(repmat(temp',3,1),240,1);
    metadata.True3D(a:b) = temp;
    
    temp = any(bsxfun(@eq,trial_vector, [21, 22, 37,38]),2);
    temp = reshape(repmat(temp',3,1),240,1);
    metadata.TrueSlow(a:b) = temp;
    
    temp = any(bsxfun(@eq,trial_vector, [25, 26, 41, 42]),2);
    temp = reshape(repmat(temp',3,1),240,1);
    metadata.TrueFast(a:b) = temp;
    
    temp = any(bsxfun(@eq,trial_vector, [22, 26, 21, 25]),2);
    temp = reshape(repmat(temp',3,1),240,1);
    metadata.TrueLeft(a:b) = temp;
    
    temp = any(bsxfun(@eq,trial_vector, [38, 42, 37, 41]),2);
    temp = reshape(repmat(temp',3,1),240,1);
    metadata.TrueRight(a:b) = temp;
    
    metadata.mriRuns(a:b) = ii;
end

%% Random balanced CVs (this is probably a bad idea because of autocorrelation).
% cv = zeros(2880/3,1);
% temp = metadata.True2D(mod(1:2880,3)==1);
% cv(temp) = mod(1:sum(temp),10); 
% 
% temp = metadata.True3D(mod(1:2880,3)==1);
% cv(temp) = mod(1:sum(temp),10); 
% 
% temp = metadata.TrueSlow(mod(1:2880,3)==1);
% cv(temp) = mod(1:sum(temp),10); 
% 
% temp = metadata.TrueFast(mod(1:2880,3)==1);
% cv(temp) = mod(1:sum(temp),10); 
% 
% cv = reshape(repmat(cv',3,1),2880,1);
% 
% CV = bsxfun(@eq, cv, (0:9));

%% CVs based on sessions (better, in that it avoids autocorrelation problem).
CV = bsxfun(@eq, metadata.mriRuns, (1:12));

%% cleanup
clear a b ii temp trial_vector vector_path


