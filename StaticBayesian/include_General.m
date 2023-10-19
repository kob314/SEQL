curren_dir  = pwd;
idcs   = strfind(curren_dir,'\');
parent_dir = curren_dir(1:idcs(end)-1);
addpath([parent_dir,'\General'])