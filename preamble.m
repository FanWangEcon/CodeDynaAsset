function [st_root_path] = preamble()
% Remote Path
st_root_path = 'C:/Users/fan/CodeDynaAsset/';
rmpath(genpath(st_root_path))
addpath(genpath(st_root_path))
end