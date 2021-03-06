%% Initialize a new mrVista session
%  The current directory will become a mrVista session. To initialize, we
%  need three things:
%   1. The inplane anatomoy (in our case, this is OutBrick_run_004_resampled.nii
%   2. The inplane functionals (OutBrick_run_{007-018}.nii, in that order.
%   3. The vAnatomy file. This is used across studies, and so might reside
%   somewhere else, physically. As a rule, there should be an Anatomy
%   directory within each mrVista session, which is a symbolic link to a
%   shared Anatomy folder.
nscans = 12;
opts = mrInitDefaultParams;
opts.inplane = 'Raw/OutBrick_run_004_resampled.nii';
opts.functionals = strsplit(sprintf('Raw/OutBrick_run_%03d.nii:',7:18),pathsep);
opts.functionals(end) = []; % the last one will always be empty, so drop it.
opts.keepFrames = repmat([12,-1],12,1); % Drop first 12 frames; keep the rest.
opts.vAnatomy = 'Anatomy/vAnatomy.dat'; % This doesn't actually set a field.

% These flags will ensure that motion correction and coherence analysis are
% done.
opts.sliceTimingCorrection = false;
opts.motionComp = 3; % Both within and between scan compensation (rigid body).
opts.applyCorAnal = true;

% These are the parameters references for the coherence analysis.
coParams = struct(...
	'blockedAnalysis', true, ...
	'detrend', true, ...
	'inhomoCorrect', true, ...
	'temporalNormalization', false, ...
	'nCycles', 80);
opts.coParams = cell(1,nscans);
for s = 1:nscans
    opts.coParams{s} = coParams;
end

% Turn off waitbars.
setpref('VISTA','verbose',0)

% Run all the things.
mrInit(opts);

% For some reason, this does not seem to be run properly by mrInit. So,
% we'll tell it a second time. Don't make me say it a third time, Mr.
% Vista.
INPLANE{1} = computeCorAnal(INPLANE{1}, opts.applyCorAnal, 1);


% Turn on waitbars.
setpref('VISTA','verbose',1)

setVAnatomyPath(opts.vAnatomy)
