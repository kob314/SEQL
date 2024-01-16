curren_dir  = pwd;
idcs   = strfind(curren_dir,'\');
parent_dir = curren_dir(1:idcs(end)-1);
addpath([parent_dir,'\General'])
addpath([parent_dir,'\numpyro_out\psychom'])

o.color = [0.4 0.55 0.25; .45 .25 .5;.9 .35 .4; .2 .3 .6; 0 0 0; .86 .45 .1; 0.93 0.74 0.35];
o.color([14 17],:) = [.6 .35 .2; [1 1 1]*.3];
