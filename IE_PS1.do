/**********************************************************************
* Homar A. Maurás Rodríguez
* Impact Evaluation – Problem Set 1
**********************************************************************/

///////////////////////////////////////////////////////////////////////////////
// QUESTION 1.1.1.d
// Simulate distribution of RED balloons and derive critical values
///////////////////////////////////////////////////////////////////////////////

cap prog drop p111d
prog def p111d, rclass
    args samplesize
    clear
    set obs 5
    gen pass = runiform() > 0.9
    return scalar N = `c(N)'
    count if pass
    return scalar pass = `r(N)'
end

simulate total = r(N) pass = r(pass), reps(10000) nodots: p111d 100
xtile pct = pass, n(1000)

gen one = 1
gen zero = 0
sort pass

tw (rarea zero one pass if pct < 25, yaxis(2) lw(none) fc(red%50)) ///
   (rarea one zero pass if pct > 975, yaxis(2) lw(none) fc(red%50)) ///
   (histogram pass, fc(gray) barw(0.9) lc(none) start(0) w(1)), ///
   yscale(off axis(2)) ytit("", axis(2)) yscale(alt) ///
   xtitle("Number Passed out of 100")

///////////////////////////////////////////////////////////////////////////////
// QUESTION 1.1.2.d
///////////////////////////////////////////////////////////////////////////////

cap prog drop p112d
prog def p112d, rclass
    args samplesize
    clear
    set obs 5
    gen pass = runiform() > 0.8
    return scalar N = `c(N)'
    count if pass
    return scalar pass = `r(N)'
end

simulate total = r(N) pass = r(pass), reps(10000) nodots: p112d 100
xtile pct = pass, n(1000)

gen one = 1
gen zero = 0
sort pass

tw (rarea zero one pass if pct < 25, yaxis(2) lw(none) fc(red%50)) ///
   (rarea one zero pass if pct > 975, yaxis(2) lw(none) fc(red%50)) ///
   (histogram pass, fc(gray) barw(0.9) lc(none) start(0) w(1)), ///
   yscale(off axis(2)) ytit("", axis(2)) yscale(alt) ///
   xtitle("Number Passed out of 100")

///////////////////////////////////////////////////////////////////////////////
// QUESTION 1.1.4 – Exact Binomial Calculation
///////////////////////////////////////////////////////////////////////////////

cap prog drop binomial_calc_exact
prog def binomial_calc_exact
    syntax, trials(integer) prob(real) successes(integer)
    local probability = binomialp(`trials', `successes', `prob')
    di "P(X = `successes') with n = `trials' and p = `prob' is: " `probability'
end

binomial_calc_exact, trials(5) prob(0.1) successes(0)

///////////////////////////////////////////////////////////////////////////////
// QUESTION 1.2 – Simulating Hypothesis Testing with Different Pass Rates
///////////////////////////////////////////////////////////////////////////////

cap prog drop p121
prog def p121, rclass
    args samplesize
    clear
    set obs `samplesize'
    gen pass = runiform() > 0.5
    return scalar N = `c(N)'
    count if pass
    return scalar pass = `r(N)'
end

simulate total = r(N) pass = r(pass), reps(10000) nodots: p121 25
xtile pct = pass, n(1000)

gen one = 1
gen zero = 0
sort pass

tw (rarea zero one pass if pct < 25, yaxis(2) lw(none) fc(red%50)) ///
   (rarea one zero pass if pct > 975, yaxis(2) lw(none) fc(red%50)) ///
   (histogram pass, fc(gray) barw(0.9) lc(none) start(0) w(1)), ///
   yscale(off axis(2)) ytit("", axis(2)) yscale(alt) ///
   xtitle("Number Passed out of 25") xla(6(2)25)

///////////////////////////////////////////////////////////////////////////////
// Repeat for Pass Rate 55%
///////////////////////////////////////////////////////////////////////////////

cap prog drop p121
prog def p121, rclass
    args samplesize
    clear
    set obs `samplesize'
    gen pass = runiform() > 0.55
    return scalar N = `c(N)'
    count if pass
    return scalar pass = `r(N)'
end

simulate total = r(N) pass = r(pass), reps(10000) nodots: p121 25
xtile pct = pass, n(1000)

gen one = 1
gen zero = 0
sort pass

tw (rarea zero one pass if pct < 25, yaxis(2) lw(none) fc(red%50)) ///
   (rarea one zero pass if pct > 975, yaxis(2) lw(none) fc(red%50)) ///
   (histogram pass, fc(gray) barw(0.9) lc(none) start(0) w(1)), ///
   yscale(off axis(2)) ytit("", axis(2)) yscale(alt) ///
   xtitle("Number Passed out of 25") xla(6(2)25)

///////////////////////////////////////////////////////////////////////////////
// QUESTION 1.2.3 – Binomial Test Over Grid of Hypotheses
///////////////////////////////////////////////////////////////////////////////

local flips = 25
local heads = 9
local alpha = 0.05

clear
set obs 99
gen p_heads = _n / 100
gen p_value = .

forvalues i = 1/99 {
    local prob = `i' / 100
    quietly bitesti `flips' `heads' `prob', d
    replace p_value = r(p) in `i'
}

gen not_rejected = p_value >= `alpha'
list p_heads p_value not_rejected if not_rejected == 1

summarize p_heads if not_rejected == 1
display "The confidence interval for the probability of heads is from " ///
    r(min)*100 "% to " r(max)*100 "%"

///////////////////////////////////////////////////////////////////////////////
// QUESTION 1.2.4 – Z-test Based on Varying p
///////////////////////////////////////////////////////////////////////////////

clear
set obs 100
gen p = (_n - 1) / 99

local x = 18
local n = 47

gen z = (`x' - `n'*p) / sqrt(`n'*p*(1 - p))
gen p_value = 2 * (1 - normal(abs(z)))

twoway (line p_value p, sort), ///
    ytitle("P-Value") ///
    xtitle("Probability of Success (p)") ///
    title("P-Value vs Probability of Success") ///
    ylabel(0(0.1)1) ///
    xlabel(0(0.1)1)

///////////////////////////////////////////////////////////////////////////////
// QUESTION 2.1 – Simulating Power & Significance in Multiple Testing
///////////////////////////////////////////////////////////////////////////////

cap prog drop p1
prog def p1, rclass
    args ptrue
    clear
    set obs 1

    gen true = runiform() < `ptrue'
    return scalar true = true[1]

    gen effect = runiform() * 0.1
    replace effect = 0 if true == 0

    expand `=ceil(runiform(100,1000))'
    gen pass = runiform() + effect > 0.01
    count if pass

    quietly bitesti `c(N)' `r(N)' 0.5, d
    return scalar significant = r(p) < 0.1
end

simulate true = r(true) significant = r(significant), reps(1000) nodots: p1 0.05

// Post-simulation inference
count if significant == 1
local p_significant = r(N) / _N

count if true == 1 & significant == 1
local p_significant_true = r(N)
count if true == 1
local total_true = r(N)
local p_significant_true = `p_significant_true' / `total_true'

count if true == 0 & significant == 0
local p_not_significant_not_true = r(N)
count if true == 0
local total_not_true = r(N)
local p_not_significant_not_true = `p_not_significant_not_true' / `total_not_true'

count if true == 1 & significant == 0
local p_not_significant_true = r(N)
local p_not_significant_true = `p_not_significant_true' / `total_true'

count if true == 0 & significant == 1
local p_significant_not_true = r(N)
local p_significant_not_true = `p_significant_not_true' / `total_not_true'

count if true == 0 & significant == 1
local p_not_true_significant = r(N)
count if significant == 1
local total_significant = r(N)
local p_not_true_significant = `p_not_true_significant' / `total_significant'

di "P(significant) = `p_significant'"
di "P(significant | true) = `p_significant_true'"
di "P(not significant | not true) = `p_not_significant_not_true'"
di "P(not significant | true) = `p_not_significant_true'"
di "P(significant | not true) = `p_significant_not_true'"
di "P(not true | significant) = `p_not_true_significant'"

///////////////////////////////////////////////////////////////////////////////
// QUESTION 3.1.1 – Simulate Mean Height from Normal Distribution
///////////////////////////////////////////////////////////////////////////////

cap prog drop mean_sample_height
prog define mean_sample_height, rclass
    args n
    clear
    set obs `n'
    gen height = rnormal(177, 7)
    summarize height, meanonly
    return scalar mean_height = r(mean)
end

mean_sample_height 10000
di r(mean_height)

///////////////////////////////////////////////////////////////////////////////
// QUESTION 3.1.2 – Sampling Distribution of the Mean
///////////////////////////////////////////////////////////////////////////////

cap prog drop mean_sample_height
prog define mean_sample_height, rclass
    args n
    clear
    set obs `n'
    gen height = rnormal(177, 7)
    summarize height
    return scalar mean_height = r(mean)
    return scalar standard_error = 7 / sqrt(42)
end

simulate mean = r(mean_height), reps(10000) nodots: mean_sample_height 42
gen id = _n
rename mean mean_height

gen lower_bound = 177 - 1.96 * (7 / sqrt(42))
gen upper_bound = 177 + 1.96 * (7 / sqrt(42))

di upper_bound
di lower_bound

hist mean_height, ///
    xline(`=lower_bound', lcolor(red) lwidth(medium) lpattern(dash)) ///
    xline(`=upper_bound', lcolor(red) lwidth(medium) lpattern(dash)) ///
    title("Histogram of Mean Heights with 95% CI")

tempfile sim_men
save `sim_men', replace

///////////////////////////////////////////////////////////////////////////////
// QUESTION 3.1.3 – Simulate France Sample Mean
///////////////////////////////////////////////////////////////////////////////

cap prog drop france_height
prog define france_height, rclass
    args n
    clear
    set obs `n'
    gen height = rnormal(179, 7)
    summarize height, meanonly
    return scalar mean_height_fr = r(mean)
end

france_height 42
di r(mean_height_fr)

tempfile sim_fr
save `sim_fr'

///////////////////////////////////////////////////////////////////////////////
// QUESTION 3.1.4 – Z-test Comparing France to Null
///////////////////////////////////////////////////////////////////////////////

france_height 42
scalar mean_france = r(mean_height_fr)
scalar null_mu = 177
scalar se = 7 / sqrt(42)
scalar z = (mean_france - null_mu) / se
scalar p_value = 2 * (1 - normal(abs(z)))

di "Z-statistic: " z
di "P-value: " p_value

///////////////////////////////////////////////////////////////////////////////
// QUESTION 3.1.5 – Power of France Hypothesis Test
///////////////////////////////////////////////////////////////////////////////

cap prog drop test_hypothesis
prog define test_hypothesis, rclass
    args n
    clear
    set obs `n'
    gen height_france = rnormal(179, 7)
    summarize height_france, meanonly
    scalar sample_mean_fr = r(mean)

    scalar null_mean = 177
    scalar se = 7 / sqrt(`n')
    scalar lower_bound = null_mean - 1.96 * se
    scalar upper_bound = null_mean + 1.96 * se

    if (sample_mean_fr >= lower_bound & sample_mean_fr <= upper_bound) {
        return scalar reject_null = 0
    }
    else {
        return scalar reject_null = 1
    }
end

simulate reject_null = r(reject_null), reps(1000) nodots: test_hypothesis 42

count if reject_null == 0
scalar fail_to_reject = r(N)

scalar power = 1 - fail_to_reject / 1000
di "Power of the test: " power

///////////////////////////////////////////////////////////////////////////////
// QUESTION 3.1.6 – India Simulation and Hypothesis Test
///////////////////////////////////////////////////////////////////////////////

cap prog drop india_height
prog define india_height, rclass
    args n
    clear
    set obs `n'
    gen height = rnormal(165, 7)
    summarize height, meanonly
    return scalar mean_height_india = r(mean)
end

india_height 42
di r(mean_height_india)

scalar mean_india = r(mean_height_india)
scalar null_mu = 177
scalar se = 7 / sqrt(42)
scalar z = (mean_india - null_mu) / se
scalar p_value = 2 * (1 - normal(abs(z)))

di "Z-statistic: " z
di "P-value: " p_value

cap prog drop test_hypothesis_in_sem
prog define test_hypothesis_in_sem, rclass
    args n
    clear
    set obs `n'
    gen height_india = rnormal(165, 7)
    summarize height_india, meanonly
    scalar sample_mean_in = r(mean)

    scalar null_mean = 177
    scalar se = 7 / sqrt(`n')
    scalar lower_bound = null_mean - 1.96 * se
    scalar upper_bound = null_mean + 1.96 * se

    if (sample_mean_in >= lower_bound & sample_mean_in <= upper_bound) {
        return scalar reject_null = 0
    }
    else {
        return scalar reject_null = 1
    }
end

simulate reject_null = r(reject_null), reps(1000) nodots: test_hypothesis_in_sem 42

count if reject_null == 0
scalar fail_to_reject = r(N)
scalar power = 1 - fail_to_reject / 1000
di "Power of the test: " power

tempfile sim_in
save `sim_in'

///////////////////////////////////////////////////////////////////////////////
// QUESTION 3.1.7 – Simulating Exaggeration Bias
///////////////////////////////////////////////////////////////////////////////

cap prog drop simulate_exaggeration
prog define simulate_exaggeration, rclass
    args n true_mean
    clear
    set obs `n'
    gen height = rnormal(`true_mean', 7)
    summarize height, meanonly
    scalar sample_mean = r(mean)
    scalar null_mean = 177
    scalar effect_size = sample_mean - null_mean
    scalar se = 7 / sqrt(`n')
    scalar z = (sample_mean - null_mean) / se
    scalar p_value = 2 * (1 - normal(abs(z)))

    return scalar effect_size = effect_size
    if (p_value < 0.05) {
        return scalar reject_null = 1
    }
    else {
        return scalar reject_null = 0
    }
end

simulate effect_size = r(effect_size) reject_null = r(reject_null), reps(1000) nodots: simulate_exaggeration 42 165
simulate effect_size = r(effect_size) reject_null = r(reject_null), reps(1000) nodots: simulate_exaggeration 42 175

gen significant_165 = reject_null == 1
gen significant_effect_size_165 = effect_size if significant_165 == 1

gen significant_175 = reject_null == 1
gen significant_effect_size_175 = effect_size if significant_175 == 1

summarize significant_effect_size_165
scalar avg_effect_size_165 = r(mean)
di "Average effect size for true mean 165: " avg_effect_size_165

summarize significant_effect_size_175
scalar avg_effect_size_175 = r(mean)
di "Average effect size for true mean 175: " avg_effect_size_175
