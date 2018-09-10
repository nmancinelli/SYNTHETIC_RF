%% Script to load receiver function mat files from racepoint
%   The code does a median stack and bootstrap resamples for error
%   estimation. -N.J. Mancinelli August 2018
%
%   Note this is quite a bit simpler than the python script NJM
%   used to generate the stacks in the GRL paper, but it gets 80%
%   of the job done.
%
%

%% File Information
clear;

path = '/Users/mancinelli/PROJECTS/ARRAY_STACK/ReceiverFunctions/CRATONS/CRATON';

%network = 'II';
%station = 'NRIL';

network = 'CN';
station = 'YKW3';

lowT  = 2;
highT = 100;
DeconMethod = 'ETMTM';
filename = sprintf('RF_Depth_%ds_%ds_%s_UsePostfilter_0.mat',lowT,highT,DeconMethod);

fullpath2file = sprintf('%s/%s/%s/%s',path,network,station,filename);

%% Load File

load(fullpath2file)

%% Bootstrap depth migrated RFs

nboot = 100;
[nrfs,npts] = size(rfs);

random_samples = zeros(nboot,npts);

for iboot = 1:nboot
    
    random_indeces = ceil(rand(1,nrfs)*nrfs);
    
    %random_indeces = [9,3,4];
    
    stack_median = nanmedian(rfs(random_indeces,:));
    
    random_samples(iboot,:) = stack_median;
    
end

bootmean = mean(random_samples)';
bootstd  = std(random_samples)';

%% Draw figure

nsigma = 2.0;

hold on

f1 = bootmean-nsigma*bootstd;
f2 = bootmean+nsigma*bootstd;

plot(f1, 0:npts-1,'black');
plot(f2, 0:npts-1,'black');

xlim([-0.2,0.2]);
ylim([0,300]);

set(gca,'YDir','reverse');

title(sprintf('%s.%s',network,station));
ylabel('Depth (km)');
xlabel('RF Amplitude (Relative to Parent)');