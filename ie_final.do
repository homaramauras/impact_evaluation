/************************************************************
Homar A. Maurás Rodríguez – Final Project
Impact Evaluation & Hypothesis Testing (Stata)
************************************************************/

clear all
set more off
set seed 1099

/**************************************************************************
1. Hypothesis Testing: Visualizing Rejection Region for Normal Distribution
**************************************************************************/

// Generate a sample from the null distribution
gen x = rnormal(5, 0.5)

// Plot theoretical density and empirical kernel density
twoway (kdensity x, color(black)) ///
       (function y = normalden((x - 5)/0.5)*2, range(3.5 6.5) color(black)), ///
       xline(4.02 5.98, lcolor(red) lpattern(dash)) ///
       xlabel(4.02 5.98, grid) ///
       title("Two-Sided Rejection Region for Hypothesis Test") ///
       ytitle("Density") xtitle("Sample Mean") legend(off)

clear

// Generate data under H0 and Ha
gen x_h0 = rnormal(5, 0.5)
gen x_ha = rnormal(6.25, 0.5)

// Plot empirical and theoretical densities for both H0 and Ha
twoway ///
    (kdensity x_h0, lcolor(black)) ///
    (kdensity x_ha, lcolor(red)) ///
    (function y = normalden((x - 5)/0.5)*2, range(3.5 6.5) lcolor(black)) ///
    (function y = normalden((x - 6.25)/0.5)*2, range(4.5 7.5) lcolor(blue)), ///
    xline(4.02 5.98, lcolor(red) lpattern(dash)) ///
    xlabel(4.02 5.98 6.25, grid) ///
    text(0.04 6.25 "Sample Mean = 6.25", place(nw) color(red)) ///
    text(0.02 5.5 "p-value = 0.0124", place(sw) color(red)) ///
    title("Rejection Region with H0 and Ha Distributions") ///
    ytitle("Density") xtitle("Sample Mean") ///
    legend(order(1 "H0 Distribution" 2 "Ha Distribution") position(5))

/*******************************************
2. Simulation: Type I Error and Power (IV)
*******************************************/

local sims 2000
local N 100
local type1_error 0
local beta_null 0
local beta_alt 0
local power 0

forval i = 1/`sims' {
    clear
    set obs `N'

    gen ability = runiform()
    gen og_degree = ability > 0.5
    gen treatment = (runiform() > 0.5) & (og_degree == 0)
    gen new_degree = og_degree | treatment
    gen error_term = rnormal(0, 300)

    // Alternative outcome
    gen wage_alt = 1000 + 500 * ability + 500 * new_degree + error_term
    ivregress 2sls wage_alt (new_degree = treatment)
    local beta_alt = `beta_alt' + _b[new_degree]
    if abs(_b[new_degree]/_se[new_degree]) > invttail(`N'-2, 0.025) {
        local power = `power' + 1
    }

    // Null outcome
    gen wage_null = 1000 + 500 * ability + error_term
    reg wage_null new_degree
    local beta_null = `beta_null' + _b[new_degree]
    if abs(_b[new_degree]/_se[new_degree]) > invttail(`N'-2, 0.025) {
        local type1_error = `type1_error' + 1
    }
}

di "Sample Size: `N'"
di "Type I Error Rate (%): " (`type1_error'/`sims')*100
di "Mean Beta under Null: " `beta_null'/`sims'
di "Mean Beta under Alternative: " `beta_alt'/`sims'
di "Bias: " (`beta_alt'/`sims') - 500
di "Power (%): " (`power'/`sims')*100

/*******************************************
3. Simulate with Program: Type I & II Errors
*******************************************/

cap prog drop q7sim
program define q7sim, rclass
    args N
    clear
    set obs `N'

    gen ability = runiform()
    gen treatment = ability < 0.5
    gen error_term = rnormal(0, 300)

    gen wage_null = 1000 + 500 * ability + error_term
    gen wage_alt = 1000 + 500 * ability + 500 * treatment + error_term

    reg wage_null treatment
    return scalar beta_null = _b[treatment]
    return scalar type1_error = (r(p) < 0.05)

    reg wage_alt treatment
    return scalar beta_alt = _b[treatment]
    return scalar bias = _b[treatment] - 500
end

simulate type1_error=r(type1_error) beta_null=r(beta_null), reps(2000) seed(1099) nodots: q7sim 1000
summarize type1_error beta_null

simulate beta_alt=r(beta_alt) bias=r(bias), reps(2000) seed(1099) nodots: q7sim 1000
summarize beta_alt bias

/*******************************************
4. Percentile-Based Treatment Assignment
*******************************************/

use "/Users/homi/Downloads/testscores.dta", clear
sort testscore0
set seed 1099

local percentiles = 100
matrix results = J(`percentiles', 2, .)

forval p = 1/`percentiles' {
    quietly sum testscore0 if _n <= `p'/`percentiles' * _N
    local cutoff = r(max)

    gen tutoring = (testscore0 <= `cutoff')
    reg testscore1 tutoring c.testscore0##c.testscore0 if !missing(testscore0, testscore1)

    matrix results[`p', 1] = _b[tutoring]
    matrix results[`p', 2] = (abs(_b[tutoring]/_se[tutoring]) > invttail(e(df_r), 0.025))
}

svmat results, names(col)
gen cutoff_percentile = _n
rename c1 beta
rename c2 reject

summ reject if !missing(reject)
local rejection_rate = (r(sum)/`percentiles')*100
di "Rejection Rate for the Null Hypothesis: `rejection_rate'%"


