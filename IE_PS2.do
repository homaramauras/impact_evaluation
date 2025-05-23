* Homar A. Maurás Rodríguez
* Impact Evaluation Problem Set

///////////////////////////////////////////////////////////////////////////////
// WARMUP
///////////////////////////////////////////////////////////////////////////////

use "/Users/homi/Downloads/regression-tools.dta", clear

* Basic regression
reg y x
return list
matlist r(table)

* Manual calculation of predicted values and residuals
reg y x
matrix b = e(b)
matrix list b

local beta0 = b[1,2]
local beta1 = b[1,1]

gen y_hat = `beta0' + `beta1' * x
gen residual = y - y_hat
gen residual_sq = residual^2

* Alternative using predict
predict xb, xb
predict resid, resid

* Compare manual and predicted results
gen y_diff = y_hat - xb
order xb y_hat y_diff resid residual
sum y_diff residual resid, detail

* Visual checks
scatter y_hat xb, mcolor(blue)
scatter resid residual, mcolor(red)
reg y xb
reg resid residual

gen r_diff = residual - resid
sum y_diff r_diff, detail

* Final plot with regression line
twoway (scatter y x, mcolor(black)) ///
       (line xb x, lcolor(red)), ///
       title("Regression Fit: y on x")

///////////////////////////////////////////////////////////////////////////////
// PROBLEM 1 – Simulating T-Distributions /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

cap prog drop sim_t_distr
prog define sim_t_distr
    args N
    clear
    set obs `N'
    gen outcome = rnormal(0,1)

    summarize outcome, detail
    local sample_mean = r(mean)
    local sample_sd = r(sd)
    local t_stat = `sample_mean' / (`sample_sd' / sqrt(`N'))

    list outcome
    di "Sample Mean: " `sample_mean'
    di "Sample SD: " `sample_sd'
    di "T-Stat: " `t_stat'
end

* Example usage
sim_t_distr 10

///////////////////////////////////////////////////////////////////////////////
// PROBLEM 1c – Repeat 10,000 Times and Plot //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

cap prog drop sim_t_stat
prog define sim_t_stat, rclass
    args N
    clear
    set obs `N'
    gen outcome = rnormal(0,1)

    summarize outcome, detail
    local sample_mean = r(mean)
    local sample_sd = r(sd)
    local t_stat = `sample_mean' / (`sample_sd' / sqrt(`N'))
    return scalar t_stat = `t_stat'
end

* Simulate and plot
simulate t_stat = r(t_stat), reps(10000) nodots: sim_t_stat 3
summarize t_stat, detail
twoway (histogram t_stat, density), title("T-Stat Distribution, N=3")

///////////////////////////////////////////////////////////////////////////////
// PROBLEM 1d – Loop Across Sample Sizes //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

clear
save "t_statistics_combined.dta", replace emptyok

foreach N in 3 5 10 30 100 1000 {
    simulate t_stat = r(t_stat), reps(10000) nodots: sim_t_stat `N'
    gen sample_size = `N'
    append using "t_statistics_combined.dta"
    save "t_statistics_combined.dta", replace
}

use "t_statistics_combined.dta", clear

twoway ///
    (kdensity t_stat if sample_size == 5, lcolor(blue)) ///
    (kdensity t_stat if sample_size == 10, lcolor(red)) ///
    (kdensity t_stat if sample_size == 30, lcolor(green)) ///
    (kdensity t_stat if sample_size == 100, lcolor(black)) ///
    (kdensity t_stat if sample_size == 1000, lcolor(magenta)), ///
    legend(order(1 "N=5" 2 "N=10" 3 "N=30" 4 "N=100" 5 "N=1000")) ///
    title("PDFs of T-Statistics by Sample Size")

///////////////////////////////////////////////////////////////////////////////
// PROBLEM 2 – Simulate Observed Wage Outcomes ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

cap prog drop p2g
prog define p2g
    clear
    set obs 3000

    gen ability_score = runiform()
    gen p_grad = .5 * ability_score
    gen error_term = rnormal(0, 10000)
    gen junior_wage = 25000 + 15000 * ability_score + error_term
    gen senior_wage = 35000 + 15000 * ability_score + error_term
    gen graduated = (runiform() < ability_score)
    gen p_senior = .2 + .25 if graduated
    gen senior_job = (runiform() < p_senior)
    gen observed_wage = senior_job * senior_wage + (1 - senior_job) * junior_wage

    sort observed_wage
    gen wage_rank = _n
    gen percentile_rank = (wage_rank - 1) / (_N - 1)

    gen marriage_prob = 0.25 + 0.15 * graduated + 0.5 * percentile_rank
    gen married = (runiform() < marriage_prob)
end

p2g
list ability_score graduated observed_wage percentile_rank marriage_prob married in 1/10

///////////////////////////////////////////////////////////////////////////////
// PROBLEM 3 – Simulate Regressions and Assess Bias ///////////////////////////
///////////////////////////////////////////////////////////////////////////////

cap prog drop p4
prog define p4, rclass
    clear
    set obs 3000

    gen ability_score = runiform()
    gen p_grad = .5 * ability_score
    gen error_term = rnormal(0, 10000)
    gen junior_wage = 25000 + 15000 * ability_score + error_term
    gen senior_wage = 35000 + 15000 * ability_score + error_term
    gen graduated = (runiform() < ability_score)
    gen p_senior = .2 + .25 if graduated
    gen senior_job = (runiform() < p_senior)
    gen observed_wage = senior_job * senior_wage + (1 - senior_job) * junior_wage

    sort observed_wage
    gen wage_rank = _n
    gen percentile_rank = (wage_rank - 1) / (_N - 1)
    gen marriage_prob = 0.25 + 0.15 * graduated + 0.5 * percentile_rank
    gen married = (runiform() < marriage_prob)

    args use_ability use_job use_married

    if `use_ability' & `use_job' & `use_married' reg observed_wage graduated ability_score senior_job married
    else if `use_ability' & `use_job' reg observed_wage graduated ability_score senior_job
    else if `use_ability' & `use_married' reg observed_wage graduated ability_score married
    else if `use_job' & `use_married' reg observed_wage graduated senior_job married
    else if `use_ability' reg observed_wage graduated ability_score
    else if `use_job' reg observed_wage graduated senior_job
    else if `use_married' reg observed_wage graduated married
    else reg observed_wage graduated

    return scalar beta_hat = _b[graduated]
end

* Run all combinations of covariates
forvalues a = 0/1 {
    forvalues j = 0/1 {
        forvalues m = 0/1 {
            simulate beta_hat = r(beta_hat), reps(10000) nodots: p4 `a' `j' `m'
            save sim_results_`a'_`j'_`m'.dta, replace
        }
    }
}

* Load summary stats for each model
use summary_stats.dta, clear

* Create a numeric ID for plotting
gen model_id = _n  

* Optional: Add value labels for readability
label define model_labels 1 "None" 2 "Married" 3 "Job" 4 "Job+Married" ///
                         5 "Ability" 6 "Ability+Married" ///
                         7 "Ability+Job" 8 "All"
label values model_id model_labels

* Plot mean beta_hat and 95% CI for each model
twoway (scatter mean_beta_hat model_id, mcolor(blue) msymbol(circle)) ///
       (rcap ci_lower ci_upper model_id, lcolor(black)) ///
       , yline(5000, lcolor(red) lpattern(dash)) ///
         title("Estimated Effect of College on Wages") ///
         ytitle("Beta-Hat Estimate") ///
         xtitle("Model Specification") ///
         xlabel(1(1)8, valuelabel angle(45)) ///
         legend(off)
