Find the best preprocessing method for NIR calibration (using Matlab)
====

Given a spetrum dataset composed of spectrum data X (n by m) and corresponding property vector Y (n by 1), you can use this scrpit to optimize the preprocessing methods. Here is an example:




>  Load data, spectra is a built-in dataset in Matlab

```Matlab
load spectra;X=NIR;Y=octane;
```

> Step 1: Optimize the parameter w,p,d for Savitzky-Golay smoothing. 
By default, you'll get a plot with 4 rows and 3 columns, where the row represent different order of fitted polynomial (0\~3).
The column represent different order of derivative (0\~2), and the window length were given in each subplot

```Matlab
sg_figure_handle = sg_enumerate(X,Y);
```

![](https://github.com/jiaoyiping630/spectrum-preprocessing/blob/master/images/SG_para_opt_for_octane.jpg) 

> From the figure you got, you probably find that,
the best para for SG smoothing is p = 2 or 3, w = 23, which got 0.22405 mean RMSEP. 
And the best para for SG1D is p = 1 or 2, w = 11; 
and the best para for SG2D is p = 2 or 3, w = 25. 

> Step 2: Optimize the parameter C (iterations) for OSC and NAS. 
You'll get two subplots, corresponding to OSC and NAS respectively.

```Matlab
osc_nas_figure_handle=osc_nas_enumerate(X,Y);
```

![](https://github.com/jiaoyiping630/spectrum-preprocessing/blob/master/images/OSC_NAS_para_opt_for_octane.jpg) 

> In the figure you got, C = 0 corresponding to using the raw spectrum. 
Probably you will find that, OSC and NAS hardly lead to performance gain in PLSR model. 
When the iterations goes higher (more 'irrelavent' signal were removed), the RMSEP tend to rise. 
Nevertheless, you can set set C = 2 to observe the interaction with other preprocessing schemes. 

> Step 3: Enumerate 108 (can be defined flexible) preprocessing schemes.
In the previous two steps, you might know that:<br>  
The best para for:<br>  
SG :    w = 23, p = 2 (either 2 or 3 doesn't matter, they are equivalent)<br>  
SG1D:   w = 11, p = 1<br>  
SG2D:   w = 25, p = 2<br>  
The practical para (C = 0 is trivial, it's equivalent to using raw spectrum ) for:<br>  
OSC:    C = 2<br>  
NAS:    C = 2<br>  
    
> Thus you define a struct as :
```Matlab
default_paras=struct;
default_paras.sg_w=11;      %   Note that we pass half width here (width = 23 -> halfwidth = 11)
default_paras.sg_p=2;
default_paras.sg1d_w=5;     %   width = 11 -> halfwidth = 5
default_paras.sg1d_p=1;
default_paras.sg2d_w=12;    %   width = 25 -> halfwidth = 12
default_paras.sg2d_p=2;
default_paras.osc_c=2;
default_paras.nas_c=2;
```

> And you can enumerate numerous preprocessing methods by:
```Matlab
numerous_figure_handle = preprocess_enumerate(X,Y,default_paras);
```

![](https://github.com/jiaoyiping630/spectrum-preprocessing/blob/master/images/Enumerate_for_octane.jpg) 

> You'll get a boxplot, where the RMSEP of each trial were given. 
An ideal method produces low RMSEP, and performs stable. In this case, the optimization of preprocessing has less performance gain. 
You could try the corn dataset (http://www.eigenvector.com/data/Corn/index.html) and you will find significant different even if only SG is considered.


Please cite this paper if you found this repositary is useful:
Jiao, Y, Li, Z, Chen, X, Fei, S. Preprocessing methods for near‐infrared spectrum calibration. Journal of Chemometrics. 2020; 34:e3306. https://doi.org/10.1002/cem.3306
