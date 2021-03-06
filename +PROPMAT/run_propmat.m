function [traces,tt,status,cmdout] = run_propmat(LAYmodel,ID,ph,samprate,inc,synthperiod,nsamps,cutf)
% [traces,tt,LAYmodel1D] = run_propmat(model,ID,ph,mindV,samprate,inc,synth period,nsamps,cutf)
% 
% Function to run the propagator matrix code for a given layerised model. 

 
if nargin < 2 || isempty(ID)
    ID = 'eg';
end
if nargin < 3 || isempty(ph)
    ph= 'Ps';
end
if nargin < 4 || isempty(samprate)
    samprate = 10;
end
if nargin < 5 || isempty(inc)
    inc = 5;
end
if nargin < 6 || isempty(synthperiod)
    synthperiod = 1;
end
if nargin < 7 || isempty(nsamps)
    nsamps = 2^12; % must be power of 2, MAKE SURE FORTRAN ARRAY IS BIG ENOUGH!
end
if nargin < 8 || isempty(cutf)
    cutf = ceil((4/synthperiod)*(nsamps/samprate));
    % cutf is not actually the frequency, it is the multiple of the
    % fundamental frequency that is the actual Nyquist. Since we want a
    % nyquist at least half of the input period (4/synthperiod), we need
    % cutf = fNyq / f_fund      where    f_fund = samprate/nsamps
end

 
if ~all(unique(factor(nsamps))==2)
    error('Nsamps must be some power of 2')
end

 
% remaining parms
obsdist = 0;
ocomps = 2; % 1 is [x,y,z], 2 is [r,t,z]

%% filenames
if ~ischar(ID), ID = num2str(ID);end
modfile = [ID,'.mod'];
execfile = [ID,'.cmd'];
odatfile = [ID,'.rtz'];
ifile = [ID,'_synth.in'];
ofile0 = [ID,'_synth.out0'];
ofile1 = [ID,'_synth.out1'];
ofile2 = [ID,'_synth.out2'];

%% =======================================================================

if strcmp(ph,'Ps')
    Vbot = LAYmodel.Vp(end);
elseif strcmp(ph,'Sp')
    Vbot = LAYmodel.Vs(end);
end
nlay = LAYmodel.nlay;

%% write to PropMatrix format
setenv('DYLD_LIBRARY_PATH', '/usr/local/bin/');
PROPMAT.writePROPMATmodfile( LAYmodel,modfile)
PROPMAT.writePROPMATparmfile(ifile, Vbot, nlay+1,nsamps,samprate,cutf) % add one layer for the halfspace
PROPMAT.writePROPMATexecfile( execfile,modfile,ifile,ofile0,ofile1,ofile2,odatfile,inc,ph,synthperiod,obsdist,ocomps)
system(['chmod +u+x ' execfile]);

%% do PropMatrix on it
[status,cmdout] = system(['./',execfile,' 1> propmat.out 2> propmat.err']);

if status ~= 0; fprintf('propmat execution error, status = %d\n', status); end


%% read PropMatrix output
[traces,tt] = PROPMAT.readPROPMATtr(odatfile);

%% delete files
% delete(execfile,odatfile,ifile,ofile1,ofile2,'synth.out');
%delete(modfile,execfile,odatfile,ifile,ofile0,ofile1,ofile2);
% % plot
% figure(2); clf, hold on
% comps = {'VERTICAL','RADIAL','TRANSVERSE'}; traces = traces(:,[3,1,2]);
% for ip = 1:3
% subplot(3,1,ip)
% plot(tt,traces(:,ip),'Linewidth',1.5)
% xlim([0 max(tt)]);
% ylabel(comps{ip},'fontsize',19,'fontweight','bold')
% end
% set(gcf,'position',[680         273        1058         825])

