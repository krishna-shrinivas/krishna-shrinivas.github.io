%   Runs 3D FISH-IF colocalization upon user-defined input
%   User needs to provide the following input:
%       o directory_name == path to images (see guide for file organization)
%       o output_dir_name == path to write output (nested within above folder)
%       o FISH_name == Name of the FISH channel (appears on plots/data)
%       o IF_name == Name of the IF channel (appears on plots/data)
%       o FISH_channel == Used to assign image file to channel 
%       o IF_channel == Used to assign image file to channel 
%
%   Examples:
%    >> call_3Danalysis_code.m --> calls params.m


tic;
addpath('source_files/');
directory_name = 'nanog_brd4/Images/';
output_dir_name = 'brd4_nanog_output';
in_params.FISH_name = 'Nanog';
in_params.IF_name = 'brd4';
in_params.FISH_channel = '561';
in_params.IF_channel = '488';
params(directory_name,output_dir_name,in_params);
time_to_complete = toc;
disp(['Simulation converged in ' num2str(time_to_complete) ' seconds']);
