function [st_root_path] = preamble(varargin)
%% PREAMBLE gets path to root and generates all path
% generate path by default, but when calling this from other functions to
% get st_root_path, do not generate path.

bl_gen_path = true;
default_params = {bl_gen_path};
[default_params{1:length(varargin)}] = varargin{:};
[bl_gen_path] = default_params{:};

%% Root Path
st_root_path = 'C:/Users/fan/CodeDynaAsset/';

%% Add to Path
if (bl_gen_path)
    rmpath(genpath(st_root_path))
    addpath(genpath(st_root_path))
end

end