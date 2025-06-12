{smcl}
{* *! version 1  29December2022}{...}
{cmd:help ebct}{right: ({browse "https://doi.org/10.1177/1536867X231196291":SJ23-3: st0726})}
{hline}

{title:Title}

{p2colset 5 13 15 2}{...}
{p2col :{cmd:ebct} {hline 2}}Estimate entropy-balancing weights for continuous
treatments, dose-response functions, and their derivatives{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 12 2}
{cmd:ebct} {help varlist:{it:covariates}} {ifin}{cmd:,} {opth treatvar(varname)}
[{opt basew(varname)} {opt samplew(varname)} {opt p(integer)}
{opt out(varname)} {opt est(string)} {opt cv} {opt bandw(real)} {opt bootstrap} {opt reps(integer)} {opt seed(integer)} {opt graph} {opt generate(stub)}]

{phang}
{it:covariates} is a list of background characteristics one wishes to control
for.


{marker description}{...}
{title:Description}

{pstd}
{opt ebct} estimates entropy-balancing weights for continuous treatments
(T{c u:}bbicke [2022]).  The routine minimizes the deviation of weights from
uniform or user-specified base weights while adhering to zero correlation and
normalization constraints.  Resulting weights can be used for nonparametric
estimation of dose-response functions (DRF) and their derivatives using
local-linear regression with an Epanechnikov kernel.

{pstd}
{cmd:ebct} provides some global summary statistics on the resulting balancing
quality from a (weighted) regression of the treatment T on the covariates X.
Reported are the regression {it:R}-squared, the overall {it:F} statistic of
the excludability of all covariates, and the corresponding {it:p}-value before
and after balancing.  For user convenience, these statistics are
stored in {cmd:e(balance)}.  Nonconvergence of the algorithm or imperfect
balance after running the {cmd:ebct} command hints toward numerical
difficulties.  If this is the case, transformation of the treatment variable
to make its distribution more symmetric, respecification of covariates, or
trimming of the sample is advised.  In addition to global balancing
indicators, the command also produces detailed balancing statistics from
linear regressions of each covariate on a fourth-order polynomial in the
treatment variable.  If covariates are truly balanced, the regression
{it:R}-squared and its corresponding {it:p}-values should be close to 0 and
1, respectively.  Results are stored in {cmd:e(balance_detail)}.  If balance
is not found to be sufficient, it is recommended that the user increase p,
that is, the order to which the treatment variable is sought to be
uncorrelated with covariates.  Lastly, {cmd:ebct} also automatically provides
summary statistics on the distribution of estimated weights.  Large weights
are likely to increase the variance of estimates.

{pstd}
The program also allows the user to use the entropy balancing for continuous
treatments weights to conveniently obtain nonparametric estimates of the DRF
or its derivatives using local linear regressions based on the Epanechnikov
kernel.  Using the {cmd:graph} option, the user plots the DRF or its
derivative against the treatment variable.  With bootstrap standard errors
available, a publication-quality graph that includes 95% confidence bands is
produced.  When available, variable labels of the outcome and the treatment
variable are used as axis titles.  Results are stored in {cmd:e()}.  In
addition, one may request the command to store estimation results as Stata
variables by using the {cmd:generate()} option.  Again, when these variables
already exist in the dataset, they will be overwritten without warning.


{marker option}{...}
{title:Options}

{phang}
{opth treatvar(varname)} specifies the continuous treatment variable.
{cmd:treatvar()} is required.

{phang}
{opt basew(varname)} specifies base weights q.  The default is
{cmd:basew(}{it:q=1/N}{cmd:)}.

{phang}
{opt samplew(varname)} specifies the sampling weights used to generate
balancing targets.

{phang}
{opt p(integer)} sets balancing targets.  {cmd:p(1)} implies that {cmd:ebct}
aims to estimate balancing weights such that the means of covariates X and
treatment variable T remain the same and that Corr(X,T)=0 after weighting.
Choosing {cmd:p(2)} additionally requires the mean of T^2 to stay the same as
well as Corr(X,T^2)=0 after weighting, and so on. The default is {cmd:p(1)}.
The maximum supported is {cmd:p(3)}.

{phang}
{opt out(varname)} specifies the name of the outcome variable for which the
DRF (or its derivative) is to be estimated.

{phang}
{opt est(string)} sets the target parameter of interest: {cmd:drf} or
{cmd:derivative}.  The default is {cmd:est(drf)}.

{phang}
{opt cv} indicates that the Epanechnikov kernel bandwidth for the estimation
of the DRF or its derivative is obtained using cross-validation.  By default,
the bandwidth is chosen to be a plugin bandwidth following Fan (1992).

{phang}
{opt bandw(real)} specifies a certain bandwidth directly instead of using
cross-validation or the plugin bandwidth.

{phang}
{opt bootstrap} indicates that standard errors shall be estimated using a
built-in bootstrap routine.

{phang}
{opt reps(integer)} sets the number of repetitions used for the bootstrap.
The default is {cmd:reps(100)}.

{phang}
{opt seed(integer)} specifies the seed used as a starting point for the random
bootstrap resampling process.  Default {cmd:seed(2222)} is used to ensure that
results are always reproducible.

{phang}
{opt graph} indicates that a graph with estimation results will be produced.
If the {cmd:bootstrap} option is specified, the graph also includes 95%
confidence intervals.  For labels of the {it:y} and {it:x} axes, the variable
labels of the outcome and the treat variable are used.  If no label exists,
variable names are used instead.

{phang}
{opt generate(stub)} allows the user to specify a {it:stub} for generating
variables that contain the estimation results.  By default, no variables are
generated.  If a {it:stub} is provided, however, several variables will be
generated depending on the options chosen.  In all cases,
{it:stub}{cmd:_weight} is generated and contains the estimated balancing
weights.  If an outcome is specified, {it:stub}{cmd:_estimates} and
{it:stub}{cmd:_grid} are generated as well.  These variables contain the
estimates of the DRF (or its derivative) and the gridpoints for which these
quantities are estimated.  If the {cmd:bootstrap} option is chosen,
{it:stub}{cmd:_stderr}, {it:stub}{cmd:_ci_low}, and {it:stub}{cmd:_ci_up} are
also generated.  These variables contain estimated standard errors and the
lower and upper bounds of the 95% confidence intervals.  If any of these
variables already exist in the dataset, they will be overwritten without
warning.


{marker examples}{...}
{title:Example using simulated data}

{pstd}
Setup{p_end}
{phang2}{cmd:. set seed 3333}{p_end}
{phang2}{cmd:. set obs 500}{p_end}
{phang2}{cmd:. generate x1 = rnormal()}{p_end}
{phang2}{cmd:. generate x2 = rnormal()}{p_end}
{phang2}{cmd:. generate treat = 0.3*(x1+x2) + rnormal(0,1)}{p_end}
{phang2}{cmd:. generate outcome = x1+x2 + rnormal(0,1)}{p_end}

{pstd}
Sample restriction{p_end}
{phang2}{cmd:. quietly summarize treat, detail}{p_end}
{phang2}{cmd:. keep if inrange(treat, r(p1), r(p99))}{p_end}

{pstd}
Estimation of DRF with bootstrap standard errors{p_end}
{phang2}{cmd:. ebct x1 x2, treatvar(treat) out(outcome) est(drf) bootstrap graph}{p_end}

{pstd}
Estimation of the derivative of the DRF with bootstrap standard errors{p_end}
{phang2}{cmd:. ebct x1 x2, treatvar(treat) out(outcome) est(derivative) bootstrap graph}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:ebct} stores the following in {cmd:e()}:

{synoptset 19 tabbed}{...}
{p2col 5 15 19 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(bandwidth)}}bandwidth used in the local linear regressions{p_end}

{p2col 5 15 19 2: Macros}{p_end}
{synopt:{cmd:e(estimand)}}either {cmd:drf} or {cmd:derivative}{p_end}

{p2col 5 15 19 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}estimated DRF or its derivative at gridpoints{p_end}
{synopt:{cmd:e(se)}}estimated standard errors of DRF or its derivative at gridpoints{p_end}
{synopt:{cmd:e(balance)}}overall balance summary statistics{p_end}
{synopt:{cmd:e(balance_detail)}}detailed balance summary statistics{p_end}
{synopt:{cmd:e(gridpoints)}}gridpoints used in the estimation{p_end}

{p2col 5 15 19 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}estimation sample{p_end}
{p2colreset}{...}


{marker references}{...}
{title:References}

{phang}
Fan, J. 1992. Design-adaptive nonparametric regression.
{it:Journal of the American Statistical Association} 87: 998–1004.
{browse "https://doi.org/10.2307/2290637":https://doi.org/10.2307/2290637}.

{phang}
T{c u:}bbicke, S. 2022. Entropy balancing for continuous treatments. 
{it:Journal of Econometric Methods} 11: 71-89.
{browse "https://doi.org/10.1515/jem-2021-0002":https://doi.org/10.1515/jem-2021-0002}.


{marker Author}{...}
{title:Author}

{pstd}Stefan T{c u:}bbicke{p_end}
{pstd}Institute for Employment Research{p_end}
{pstd}Nuremberg, Germany{p_end}
{pstd}stefan.tuebbicke@iab.de{p_end}


{title:Also see}

{p 4 14 2}
Article:  {it:Stata Journal}, volume 23, number 3: {browse "https://doi.org/10.1177/1536867X231196291":st0726}{p_end}

{p 7 14 2}
Help:  {helpb drf}, {helpb qcte} (if installed){p_end}
