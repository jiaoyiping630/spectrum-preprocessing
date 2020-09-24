%%  
%   Given a spetrum dataset composed of spectrum data X (n by m) and 
%   corresponding property vector Y (n by 1), you can optimize the preprocessing methods by following steps


%%  Load data
load spectra;X=NIR;Y=octane;    %   spectra is a built-in dataset in Matlab

%%  Step 1: Optimize the parameter w,p,d for Savitzky-Golay smoothing
%   by default, you'll get a plot with 4 rows and 3 columns
%   where the row represent different order of fitted polynomial (0~3)
%   the column represent different order of derivative (0~2)
%   the window length were given in each subplot
sg_figure_handle = sg_enumerate(X,Y);

%   From the figure you got, you probably find that,
%   the best para for SG smoothing is p = 2 or 3, w = 23, which got 0.22405 mean RMSEP.
%   And the best para for SG1D is p = 1 or 2, w = 11;
%   and the best para for SG2D is p = 2 or 3, w = 25.

%%  Step 2: Optimize the parameter C (iterations) for OSC and NAS
%   you'll get two subplots, corresponding to OSC and NAS respectively
osc_nas_figure_handle=osc_nas_enumerate(X,Y);

%   In the figure you got, C = 0 corresponding to using the raw spectrum
%   Probably you will find that, OSC and NAS hardly lead to performance gain in PLSR model.
%   When the iterations goes higher (more 'irrelavent' signal were removed), the RMSEP tend to rise.
%   Nevertheless, you can set set C = 2 to observe the interaction with other preprocessing schemes.

%%  Step 3: Enumerate 108 (can be defined flexible) preprocessing schemes
%   In the previous two steps, you might know that:
%   The best para for:
%       SG :    w = 23, p = 2 (either 2 or 3 doesn't matter, they are equivalent)
%       SG1D:   w = 11, p = 1
%       SG2D:   w = 25, p = 2
%   The practical para (C = 0 is trivial, it's equivalent to using raw spectrum ) for:
%       OSC:    C = 2
%       NAS:    C = 2
%   Thus you define a struct as :
default_paras=struct;
default_paras.sg_w=11;      %   Note that we pass half width here (width = 23 -> halfwidth = 11)
default_paras.sg_p=2;
default_paras.sg1d_w=5;     %   width = 11 -> halfwidth = 5
default_paras.sg1d_p=1;
default_paras.sg2d_w=12;    %   width = 25 -> halfwidth = 12
default_paras.sg2d_p=2;
default_paras.osc_c=2;
default_paras.nas_c=2;
%   And you can enumerate numerous preprocessing methods by:
numerous_figure_handle = preprocess_enumerate(X,Y,default_paras);

%   You'll get a boxplot like this, where the RMSEP of each trial were given.
%   An ideal method produces low RMSEP, and performs stable.

