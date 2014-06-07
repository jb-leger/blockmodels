
our %config;

$config{'MEMBERSHIP_ARG'} = "membership_type";
$config{'MODEL_ARG'} = "verbosity=6, autosave='', plotting=character(0), exploration_factor=1.5, explore_min=4, explore_max=Inf, ncores=detectCores()";
$config{'SCALAR_MODEL_ARG'} = "adj";

$config{'MEMBERSHIP_ARG_ITEM'} = 
    "\\item{membership_type}{The type of node membership, i.e. 'SBM', 'SBM_sym' or 'LBM'}";

$config{'MODEL_ARG_ITEM'} =
    "\\item{verbosity}{The verbosity level, 0 means quiet. Level 1 display the phase of reinitialization. Level 2 display the level 1 and the ascending and descending phase for the number of groups. Level 3 display the level 2 and the number current number of groups which is estimated. Level 4 display the level 3 and the steps inside the estimation. Level 5 display the level 4, the current status of parallel running jobs and the current sub-step. Level 6 display level 5 and informations about ICL criteria found. Default is level 6. This parameter can be changed by accessing to the field \$verbosity of the object.}\n".
    "\\item{autosave}{If \\var{autosave} != '', after each estimation, the model object is writed into file \\var{autosave}. The model object is readable by the function \\var{readRDS}. Use-it for long computation to allow restarting the estimation on system crash. You can use it to alanyze the partial results when the estimation is running. This parameter can be changed by accessing to the field \$autosave of the object.}\n".
    "\\item{plotting}{Control plot of ICL values while the estimation is running. If plotting==character(0) (the default), plots are done on screen, if plotting=='', no plot are done, if plotting is a filename, plots are done in this filename. This parameter can be changed by accessing the field \$plotting of the object.}".
    "\\item{exploration_factor}{Control the exploration of the number of groups. The exploration is stop when the number of groups reach exploration factor times the current maximum. By default 1.5. This parameter can be changed by accessing the field \$exploration_factor of the object.}\n".
    "\\item{explore_min}{Explore to the explore_min number of groups even if the exploration_factor rule is satisfied. By default 4. This parameter can be changed by accessing the field \$explore_min of the object.}\n".
    "\\item{explore_max}{Stop exploration after explore_max number of group in any case. By default Inf. This parameter can be changed by accessing the field \$explore_max of the object.}\n".
    "\\item{ncores}{Number of parallel jobs to launch different EM intializations. By default detectCores(). This parameter can be changed by accessing the field \$ncores of the object. This parameters is used only on Linux. Parallism is disabled on other plateform. (Not working on Windows, not tested on Mac OS, not tested on *BSD.)}\n";

$config{'SCALAR_MODEL_ARG_ITEM'} =
    "\\item{adj}{The adjacency matrix}";


$config{'COVARIATES_ARG'} = 'covariates';

$config{'COVARIATES_ARG_ITEM'} =
    "\\item{covariates}{Covariates matrix, or list of covariates matrices. Covariates matrix must have the same size than the adjacency matrix.}";
