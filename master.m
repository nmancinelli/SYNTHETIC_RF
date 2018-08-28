%% Quick-and-dirty synthetic receiver functions
%  Nicholas J. Mancinelli -- August 23, 2018
%
%  This script calls several subroutines to get the job done,
%   including a propagator matrix code that must be compiled
%   on your machine.  Look in PropMatrix/src/ for the Fortran
%   source code.  Everything else is in Matlab, thank goodness.
%
%
% 

%% Step 1: Load in velocity model from card file
%
VM = VelocityModel('STW105.txt');

%Note: For the sake of this example I've used model STW105.
%Under the hood I have extended the crust to the surface of the model.

%% Step 2: Run propmat

% Assign ray parameter and dominant period for the synthetics
ray_parameter = 0.10;
period = 7.0;

SRF = SyntheticRF(VM, ray_parameter, period);
SRF = SRF.compute_synthetic_wvfrms();

%% Step 3: Rotate the synthetics, deconvolve, and migrate to depth
SRF = SRF.rotate();
SRF = SRF.deconvolve();
SRF = SRF.migrate_to_depth();

%Note: One could add fancier windowing and masking features, if needed.

%% Step 4: Plot output
plot(SRF.rf_depth, SRF.depth,'linewidth',2)
set(gca,'Ydir','reverse')
ylabel('Depth (km)')
xlabel('Receiver Function Amplitude (Fraction of Parent Phase Amplitude)')

%% 
%saveas(gcf,'example.png')

%Note: The object SRF contains the raw synthetic waveforms, as well as the
% receiver functions in both time and depth.
