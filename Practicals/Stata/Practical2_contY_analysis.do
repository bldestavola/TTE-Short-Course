/*
---------------------------------------------------
TTE Short Course

Stata commands for Practical 2  
---------------------------------------------------
*/

*change directory to where ou ahve saved the data.
*For example
cd "S:\ICH_PPP_CENB_CEROBS\Short Course\data"

*keep a record of the results using a log file:
cap log close
log using Practical2_simul_analysis.log, replace
	
	
/*check you have read the data ok*/

*---wide format: N=10,000	
use Practical2_contY_wide,replace
des

*---long format: N=	20,000
/*
there are twice as many records in long format because there are 2 records oer individuals, corresponding to times 0 and 1
*/
use Practical2_contY_long,replace
describe
ta t


**Q1
*the arrow from L_0 to A_0 would have to be removed.


**Q2
* to summarise the data it is handier to read the data in long format because we want to count each individual only once.
use Practical2_contY_wide,replace
su
* there are 10,000 individuals

**Q3
su Y
hist Y, scheme(s1mono)
graph export Y_hist.pdf,replace
*The mean of Y is 4.04 and its SD is 1.9. The range is from -2.6 to 10.5.
* the histogram shows that it is fairly symmetrically distributioned

**Q4
ta A0
ta A0 A1, row
*Just over half (52.8\%) of the patients were prescribed treatment at time 0. Of these only over half (56.8\% were still on treatment at time 1. Of the 47.2\% who had not started treatment at time 0, 44.6\% started at time 1. Hence there is a substantial amount of treatment switching


**Q5
reg Y A0 L0  
*The estimated regression coefficient for $A_0$, controlling for $L_0$, is 1.11 (95\% CI: 1.05, 1.17) indicating a 1.1 increase in $Y$ when exposed to $A$ at time 0, condition on baseline confounder $L_0$.


**Q6
*(a) MSM-

*(b) estimate the ITT using unstandardized IPW: need Pr(A_0|L_0)
logit A0 L0 
predict den1
replace den1=1-den1 if A0==0
gen wt_A0_unst=1/den1
*ITT with IPW
su wt_A0_unst
*their mean should be 2 because we are doubling up the population
reg Y A0  [pw=wt_A0_unst]
*the estimated ITT is   1.108 (1.038,1.178), which is close to the true ITT=1.10

*(c) estimate the ITT using using standardized IPW: also need Pr(A_0)
logit A0
predict num1
replace num1=1-num1 if A0==0
gen wt_A0_stand=num1/den1
su wt_A0_stand
*their mean should be 1 because now we are standardizing  the original weights to the original population size
reg Y A0  [pw=wt_A0_stand]
*the estimated ITT is   1.108 (1.038, 1.178)


*Q7
*ITT with gcomp (using teffects)
teffects ra (Y L0) (A0)
*The same estimate is found using g-computation, although the confidence interval is tighter: 1.108 (1.050, 1.165). This is because more parametric assumptions are made by g-computation and, if appropriate,  they lead to greater precision.
 
*the two POs can be found with the option <pomeans>:
teffects ra (Y L0) (A0), pomeans


*Q8 -need to model adherence
*generate adherence indicator
gen adh=A0==A1
ta adh 
*adherence probabilities as function of A0, and L at previosu time points
logit adh A0 L0 L1 
predict ps_adh_den
replace ps_adh_den=0 if adh==0

gen wt_adh=(1/ps_adh_den)

*combining the selection into treatment and adherence weights
gen wt_adh_conf=wt_adh*wt_A0_stand
su wt_ad*
* the mean adherence weight is 1.8 which is the inverse of the prob(adherence) overall (=1/0.5611) 
*PP with IPW
reg Y A0  [pw=wt_adh_conf]
* the estimated PP is 1.986 (1.895, 2.076) which is close to the true PP=2.

 
exit
