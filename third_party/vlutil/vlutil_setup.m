function vlutil_setup
% VLUTIL_SETUP Setup path for VLUTIL toolbox
vlroot = fileparts(which('vlutil_setup')) ;
addpath(vlroot,'-begin') ;
addpath(genpath([vlroot '/toolbox']),'-begin') ;
addpath([vlroot '/mex'],'-begin') ;
fprintf('VLUTIL path loaded\n') ;
