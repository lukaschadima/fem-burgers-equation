% setup_paths.m
%
% Run this once at the start of every MATLAB session before calling
% anything in this repository. It adds all source and analysis
% subfolders to the MATLAB path so functions can find each other.
%
% Usage:
%   >> setup_paths

this_dir = fileparts(mfilename('fullpath'));

addpath(genpath(fullfile(this_dir, 'src')));
addpath(genpath(fullfile(this_dir, 'analysis')));

fprintf('Paths added. You can now run, e.g.:\n');
fprintf('  test_spatial_convergence\n');
fprintf('  test_stability_new\n');
