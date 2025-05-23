* Homar A. Maurás Rodríguez
* PPOL 6818 - Assignment Stata 4
* Simulations: Sampling Noise, Power, and Bias Adjustment

if c(username)=="homi" {
	global wd "/Users/homi/GitHub/ppol6818/ppol6818-ham86/assignments/stata_3"
}

set more off
set seed 1099

	/***********************************************
	Part 1: Power Analysis (Individual Randomization)
	***********************************************/

clear all
set obs 10000

* Outcome variable: Y ~ Normal(0,1)
gen Y = rnormal(0,1)

* 50% assigned to treatment group
gen treat = (runiform() < 0.5)

* Apply fixed average treatment effect of 0.1
replace Y = Y + 0.1 if treat == 1

* Power calculation: two-sample test with equal group sizes
power twomeans 0 0.1, sd(1) power(0.8)
scalar initial_n = r(N)

* Adjust for expected 15% attrition
scalar adjusted_n = initial_n / (1 - 0.15)
di "Adjusted sample size (accounting for attrition): " adjusted_n

* Scenario: Only 30% receive treatment -> update sample size based on imbalance
power twomeans 0 0.1, sd(1) power(0.8) nratio(0.3/0.7)

	/*******************************************
	Part 2: Power Under Cluster-Level Assignment
	*******************************************/

clear all
local cluster_size = 15
local n_clusters = 200
local rho = 0.3

* Simulate school-level effects
set obs `n_clusters'
gen school_id = _n
gen u_school = rnormal(0, sqrt(`rho'))

expand `cluster_size'
bysort school_id: gen student_id = _n
gen e = rnormal(0, sqrt(1 - `rho'))
gen Y = u_school + e

* Assign treatment to entire schools
preserve
bysort school_id (student_id): keep if _n == 1
gen treat = (runiform() < 0.5)
gen te = runiform(0.15, 0.25) if treat
replace te = 0 if missing(te)
keep school_id treat te
tempfile treatinfo
save `treatinfo'
restore

merge m:1 school_id using `treatinfo', nogen
replace Y = Y + te

* Simulate power at varying cluster sizes (m)
clear all
local sims = 500
tempname results
postfile `results' m power using "$wd/cluster_power.dta", replace

foreach m in 2 4 8 16 32 64 128 256 512 1024 {
	local rejections = 0
	forvalues s = 1/`sims' {
		clear
		set obs `n_clusters'
		gen school_id = _n
		gen u_school = rnormal(0, sqrt(`rho'))
		expand `m'
		bysort school_id: gen student_id = _n
		gen e = rnormal(0, sqrt(1 - `rho'))
		gen Y = u_school + e

		preserve
		bysort school_id (student_id): keep if _n == 1
		gen treat = (runiform() < 0.5)
		tempfile tr
		save `tr'
		restore
		merge m:1 school_id using `tr', nogen

		gen te = runiform(0.15, 0.25) if treat
		replace te = 0 if missing(te)
		replace Y = Y + te

		regress Y treat, cluster(school_id)
		if (abs(_b[treat]/_se[treat]) > invttail(e(df_r), 0.025)) {
			local rejections = `rejections' + 1
		}
	}
	local pwr = `rejections'/`sims'
	post `results' (`m') (`pwr')
}
postclose `results'

use "$wd/cluster_power.dta", clear

* Visualize power curve by cluster size
twoway line power m, ///
	title("Statistical Power by Cluster Size") ///
	xtitle("Students per School") ytitle("Power") ///
	yline(0.8, lcolor(red) lpattern(dash))

	/*************************************
	Part 3: Adjusting for Confounding Bias
	*************************************/

clear all
local sims = 500
tempname results
postfile `results' N str20 model beta using "$wd/model_betasim_output.dta", replace

foreach N in 100 200 500 1000 2000 {
	forvalues s = 1/`sims' {
		clear
		set obs `N'

		* Simulate strata groups and covariates
		gen strata = ceil(runiform()*5)
		gen x1 = rnormal() // confounder
		gen x2 = rnormal() // outcome-only predictor
		gen x3 = rnormal() // treatment-only predictor
		gen noise = rnormal()
		gen pscore = invlogit(0.5*x1 + 0.5*x3)
		gen treat = (runiform() < pscore)
		gen Y = 1*treat + 0.5*x1 + 0.3*x2 + 0.2*strata + noise

		* Model A: No covariates
		regress Y treat
		post `results' (`N') ("Model_A") (_b[treat])

		* Model B: Add confounder (x1)
		regress Y treat x1
		post `results' (`N') ("Model_B") (_b[treat])

		* Model C: Add x1 + outcome-only (x2)
		regress Y treat x1 x2
		post `results' (`N') ("Model_C") (_b[treat])

		* Model D: All three covariates
		regress Y treat x1 x2 x3
		post `results' (`N') ("Model_D") (_b[treat])

		* Model E: Fixed effects for strata
		regress Y treat i.strata
		post `results' (`N') ("Model_E") (_b[treat])
	}
}
postclose `results'

use "$wd/model_betasim_output.dta", clear
collapse (mean) est_mean=beta (sd) est_sd=beta, by(N model)

* Plot mean of beta estimates
twoway (line est_mean N if model=="Model_A", lpattern(dash)) ///
       (line est_mean N if model=="Model_B") ///
       (line est_mean N if model=="Model_C") ///
       (line est_mean N if model=="Model_D") ///
       (line est_mean N if model=="Model_E"), ///
       title("Mean Estimated Beta by Model Specification") ///
       xtitle("Sample Size (N)") ytitle("Average Estimate") ///
       yline(1, lcolor(red) lpattern(dot)) ///
       legend(order(1 "None" 2 "Confounder" 3 "x1 + x2" 4 "All Covariates" 5 "Strata FE")) ///
       name(mean_plot, replace)
graph export "$wd/plots/mean_beta_plot.png", replace

* Plot stdev of beta estimates
twoway (line est_sd N if model=="Model_A", lpattern(dash)) ///
       (line est_sd N if model=="Model_B") ///
       (line est_sd N if model=="Model_C") ///
       (line est_sd N if model=="Model_D") ///
       (line est_sd N if model=="Model_E"), ///
       title("Standard Deviation of Beta by Model Specification") ///
       xtitle("Sample Size (N)") ytitle("Standard Deviation") ///
       legend(order(1 "None" 2 "Confounder" 3 "x1 + x2" 4 "All Covariates" 5 "Strata FE")) ///
       name(sd_plot, replace)
graph export "$wd/plots/sd_beta_plot.png", replace
