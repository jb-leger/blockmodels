
our %config;

$config{'MEMBERSHIP_ARG'} = "membership_type";
$config{'MODEL_ARG'} = "verbosity=4, autosave='', plotting=character(0)";
$config{'SCALAR_MODEL_ARG'} = "adj";

$config{'MEMBERSHIP_ARG_ITEM'} = 
    "\\item{membership_type}{The type of node membership, i.e. 'SBM', 'SBM_sym' or 'LBM'}";

$config{'MODEL_ARG_ITEM'} =
    "\\item{verbosity}{The verbosity level, 0 means quiet. Level 1 display the phase of reinitialization. Level 2 display the level 1 and the ascending and descending phase for the number of groups. Level 3 display the level 2 and the number current number of groups which is estimated. Level 4 display the level 3 and the steps inside the estimation. Level 5 display the level 4 and unformations about ICL founds, and current status of parallel running jobs. Default is level 4.}\n".
    "\\item{autosave}{If \\var{autosave} != '', after each estimation, the model object is writed into file \\var{autosave}. The model object is readable by the function \\var{readRDS}. Use-it for long computation to allow restarting the estimation on system crash. You can use it to alanyze the partial results when the estimation is running.}\n".
    "\\item{plotting}{Control plot of ICL values while the estimation is running. If plotting==character(0) (the default), plots are done on screen, if plotting=='', no plot are done, if plotting is a filename, plots are done in this filename.}";

$config{'SCALAR_MODEL_ARG_ITEM'} =
    "\\item{adj}{The adjacency matrix}";


$config{'COVARIATES_ARG'} = 'covariates';

$config{'COVARIATES_ARG_ITEM'} =
    "\\item{covariates}{Covariates matrix, or list of covariates matrices. Covariates matrix must have the same size than the adjacency matrix.}";
