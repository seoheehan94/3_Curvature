%moveFilesComputeCanada.m

cd '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/stimuli';

% command = 'scp -r curvfilter/ hanseohe@beluga4.computecanada.ca:/home/hanseohe/scratch/stimuli';
% system(command);



command = 'rsync -av curvfilter/ hanseohe@beluga4.computecanada.ca:/home/hanseohe/scratch/stimuli/curvfilter';
system(command);

