% How to create a simple Taylor diagram
%
% A first example of how to create a simple Taylor diagram given one set of
% reference observations and multiple model predictions for the quantity.
% The Matlab code is kept to a minimum.
%
% This example shows how to calculate the required statistics and produce
% the Taylor diagram. All functions in the Skill Metrics Toolbox are
% designed to only work with one-dimensional arrays, e.g. time series of
% observations at a selected location. The one-dimensional data are read in
% as data structures via a mat file. The latter are stored in data
% structures in the format: ref.data, pred1.data, pred2.dat, and
% pred3.dat. The plot is written to a file in Portable Network Graphics
% (PNG) format.
%
% The reference data used in this example are cell concentrations of a
% phytoplankton collected from cruise surveys at selected locations and 
% time. The model predictions are from three different simulations that
% have been space-time interpolated to the location and time of the sample
% collection. Details on the contents of the data structures (once loaded)
% can be obtained by simply entering the data structure variable name at 
% the command prompt, e.g. 
% >> ref
% ref = 
%          data: [57x1 double]
%          date: {57x1 cell}
%         depth: [57x1 double]
%      latitude: [57x1 double]
%     longitude: [57x1 double]
%       station: [57x1 double]
%          time: {57x1 cell}
%         units: 'cell/L'
%          jday: [57x1 double]

% Author: Peter A. Rochford
%         CSS-Dynamac (Contractor)
%         NOAA/NOS/NCCOS/CCMA/COAST
%         peter.rochford@noaa.gov

% Close any previously open graphics windows
close all;

% Read in data from a mat file
load('taylor_data.mat');

% Calculate statistics for Taylor diagram
% The first array element corresponds to the reference series for the
% while the second is that for the predicted series.
taylor_stats1 = taylor_statistics(pred1,ref,'data');
taylor_stats2 = taylor_statistics(pred2,ref,'data');
taylor_stats3 = taylor_statistics(pred3,ref,'data');

% Store statistics in arrays
sdev = [taylor_stats1.sdev(1); taylor_stats1.sdev(2); ...
    taylor_stats2.sdev(2); taylor_stats3.sdev(2)];
crmsd = [taylor_stats1.crmsd(1); taylor_stats1.crmsd(2); ...
    taylor_stats2.crmsd(2); taylor_stats3.crmsd(2)];
ccoef = [taylor_stats1.ccoef(1); taylor_stats1.ccoef(2); ...
    taylor_stats2.ccoef(2); taylor_stats3.ccoef(2)];

% Produce the Taylor diagram.
%
% Note that the first index corresponds to the reference series for the
% diagram. For example sdev(1) is the standard deviation of the reference
% series and sdev(2:4) are the standard deviations of the other series.
% The value of sdev(1) is used to define the origin of the RMSD contours.
% The other values are used to plot the points (total of 3) that appear in
% the diagram.
[hp, ht, axl] = taylor_diagram(sdev,crmsd,ccoef);

% Write plot to file
writepng(gcf,'taylor1.png');
