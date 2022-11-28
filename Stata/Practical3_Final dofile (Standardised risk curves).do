global dataset = "C:\Users\micha\Desktop\Presentations\NASH (2022)\Final\Practical\Final\Final\Final"
global results = "C:\Users\micha\Desktop\Presentations\NASH (2022)\Final\Practical\Final\Final\Final\Results"

use "$dataset\Dataset (before risk curves).dta" , replace

* IP weighted risk curves 

* We run our model using functions for time (i.e. time time^2 and time^3) with interactions with the exposure
gen     IPW_A_tr=IPW_A
replace IPW_A_tr=25 if IPW_A_tr>25


logit Y_t1 A_bas male age L_bas c.time c.time#c.time c.time#c.time#c.time  A_bas#c.time A_bas#c.time#c.time A_bas#c.time#c.time#c.time  [pw=IPW_A_tr], cluster(id)



tab time
* We drop everything and we create 5 time points from scratch for each of the 5K individuals, i.e. 25K observations
drop if time ~= 0
expand 5 if time ==0 
bysort id: replace time = _n - 1		 
		
*Create 2 copies of each subject, and set outcome to missing and treatment -- use only the newobs*
expand 2 , generate(intervention) 
replace A_bas = intervention	

* Generate predicted event probabilities for each person each time point (that corresponds to 1 year) in copies*		 
* We also calculate the probability in remaining disease free, i.e. p_healthy_k
predict p_event_k, pr
gen p_healthy_k = 1-p_event_k
keep id time A_bas intervention p_healthy_k 

*Within copies, generate predicted probabilities in remaining disease free (i.e. healthy) over time*
* VERY IMPORTNANT: Remaining healthy over time is the product of conditional "remaining healthy" probabilities in each intervention*	
sort id intervention time
gen _t = time + 1
gen p_healthy = p_healthy_k if _t ==1 		
bysort id intervention: replace p_healthy = p_healthy_k*p_healthy[_t-1] if _t >1 

* Now we can calculate the probability of being at risk at each time point cumulatively (%)
gen p_risk=(1-p_healthy)*100

*Display 5-year standardized adjusted risk curves, under interventions*
*Note: since time starts at 0, time 4 is for the 5-year risk*
by intervention, sort: sum p_risk if time==4



quietly summarize p_risk if(intervention==0 & time == 4)
matrix input observe = (0,`r(mean)')
quietly summarize p_risk if(intervention==1 & time == 4)
matrix observe = (observe \1,`r(mean)')
matrix observe = (observe \3, observe[2,2]-observe[1,2]) 
matrix list observe


*Graph of standardized risk curves over time, under interventions*
*Note: since our outcome model has no covariates, we can plot p_risk directly
expand 2 if time ==0, generate(newtime)
replace p_risk  = 0 if newtime == 1
gen correct_time = 0 if newtime ==1
replace correct_time = time + 1 if newtime == 0

sort id intervention correct_time

forvalues i=0(1)5{
sum  p_risk if intervention==0 & correct_time==`i'
scalar i0_`i'=r(mean)
sum  p_risk if intervention==1 & correct_time==`i'
scalar i1_`i'=r(mean)
}

* estimate the average risk of each intervention at each time point and put it in one id, here id==1
gen     average_p_risk0=.
gen     average_p_risk1=.

forvalues i=0(1)5{
sum  p_risk if intervention==0 & correct_time==`i'
scalar i0_`i'=r(mean)
sum  p_risk if intervention==1 & correct_time==`i'
scalar i1_`i'=r(mean)
}

forvalues i=0(1)5{
local j=`i'+1
local k=`i'+7
replace average_p_risk0=i0_`i' in `j'
replace average_p_risk1=i1_`i' in `k'
}


keep if id==1

twoway (line average_p_risk0 correct_time, sort) (line average_p_risk1 correct_time, sort) if intervention > -1, ylabel(0(5)40) xlabel(0(1)5) title("Standardised Risk Curves") ytitle("Risk (%)") xtitle("Years of follow-up") /// 
       legend(label(1 "A=0") label(2 "A=1")) note("Risk A0 = 34.2% (29.7% , 38.4%)" "Risk A1 = 32.3% (30.0% , 34.8%)" "Risk difference (r(A1)-r(A0)=-1.9% (-6.5% , 3.0%)", ring(0) pos(1))


* For the estimation of the 95% confidence intervals, we used non-parametric bootstrap (minimum n of bootstraps you should use is 200, usually 500 or more is required)

capture program drop bootipw_risk_curves
program define bootipw_risk_curves , rclass
        use "$dataset\Practical TTE (before bootstrap).dta"  , clear
		preserve
		
		bsample, cluster(id) idcluster(new_id)
		
		sort new_id time
		
		/*
		the following command is very important
		
		bsample, cluster(id) idcluster(new_id)
		
		otherwise the boostrap would select random observations from each id
		
		id before bootstrap
		
		id
		1
		2
		3
		4
		5
		
		id, new_id after bootstrap
		
		seqn   newseqn
		1      1
		1      2
		1      3
		2      4
		5      5
		*/
		
		* we calculate again the the weights for the per-protocol analysis 
		
		gen NA_A0=0 if A==0 & time==0
		gen NA_A1=0 if A==1 & time==0

		forvalues i=1(1)4{

		replace NA_A0=0 if NA_A0[_n-1]==0 & A==0 & time==`i'
		replace NA_A1=0 if NA_A1[_n-1]==0 & A==1 & time==`i'

		}

		by new_id: replace NA_A0=1 if NA_A0==. & NA_A0[_n-1]==0
		by new_id: replace NA_A1=1 if NA_A1==. & NA_A1[_n-1]==0


		* we generate a variable for A on the previous time point (i.e. previous year)

		by new_id: gen A_pr=A[_n-1]



		logistic A male age L_bas L time c.time#c.time c.time#c.time#c.time if A_pr==0 & time>0 & NA_A0~=.

		predict p_NA_A0
		replace p_NA_A0=. if NA_A0~=0

		* probability of non-adherence is set to zero at baseline
		replace p_NA_A0=0 if time==0 & A==0

		* We now calculate the probability of adherence to no treatment

		gen p_A_A0= 1-p_NA_A0


		* We calculate the adherence unstabilised weights for those who initiated treatment
		* We first calculate the probability of non-adherence for treatment initiators


		logistic A male age L_bas L time c.time#c.time c.time#c.time#c.time  if A_pr==1 & time>0 & NA_A1~=.

		predict p_A_A1
		replace p_A_A1=. if NA_A1~=0

		* probability of adherence is set to one at baseline
		replace p_A_A1=1 if time==0 & A==1


		* We drop those who did not adhere (deviated treatment strategy)

		drop if NA_A0==1
		drop if NA_A1==1
		drop if NA_A0==. & NA_A1==.



		* We now calculate the IPW for adherence

		* First at each time point
		gen IPW_A0=1/p_A_A0
		gen IPW_A1=1/p_A_A1

		* We now put these weight in one variable

		gen     IPW_A=IPW_A0 if A_bas==0
		replace IPW_A=IPW_A1 if A_bas==1


		* We now multiply the probabilities at all time points

		by new_id: replace IPW_A=IPW_A*IPW_A[_n-1] if _n~=1

		* We now calculate the final (unstabilised) weights 
		* by multiplying the adherence weights IPW_A and the weights for randomisation at baseline IPW

		sum IPW_A, d
		gen IPW_A_95=r(p95)
		gen IPW_A_99=r(p99)
		gen IPW_A_max=r(max)

		* Few of our weights are rather large (IPW_final>25), however this is <1% of them (79 out of 12872)
		* For this reason, we may truncate the few exrteme weights to 25 

		* The number 25 has been chosen as a rule of thumb
		* Other researchers use as a cutoff for extreme weights from 15 to 50

		gen     IPW_A_tr=IPW_A
		replace IPW_A_tr=25 if IPW_A>25


		logit Y_t1 A_bas male age L_bas c.time c.time#c.time c.time#c.time#c.time  A_bas#c.time A_bas#c.time#c.time A_bas#c.time#c.time#c.time  [pw=IPW_A_tr], cluster(new_id)

		* We drop everything and we create 5 time points from scratch for each of the 5K individuals, i.e. 25K observations
		drop if time ~= 0
		expand 5 if time ==0 
		bysort new_id: replace time = _n - 1		 
		
		*Create 2 copies of each subject, and set outcome to missing and treatment -- use only the newobs*
		expand 2 , generate(intervention_b) 
		replace A_bas = intervention_b	

		* Generate predicted event probabilities for each person each time point (that corresponds to 1 year) in copies*		 
		* We also calculate the probability in remaining disease free, i.e. p_healthy_k
		predict p_event_k, pr
		gen p_healthy_k = 1-p_event_k
		keep new_id time A_bas intervention_b p_healthy_k IPW_A_95 IPW_A_99 IPW_A_max

		*Within copies, generate predicted probabilities in remaining disease free (i.e. healthy) over time*
		* VERY IMPORTNANT: Remaining healthy over time is the product of conditional "remaining healthy" probabilities in each interventional*	
		sort new_id intervention_b time
		gen _t = time + 1
		gen p_healthy = p_healthy_k if _t ==1 		
		bysort new_id intervention_b: replace p_healthy = p_healthy_k*p_healthy[_t-1] if _t >1 

		* Now we can calculate the probability of being at risk at each time point cumulatively (%)
		gen p_risk=(1-p_healthy)*100

		drop if time ~= 4
		bysort intervention_b: egen mean_risk_b = mean(p_risk)
		keep new_id A_bas intervention_b mean_risk_b IPW_A_95 IPW_A_99 IPW_A_max
        sort new_id intervention_b
		keep if _n<=2   /// we only need the mean values that are stored in one individual 
	
		drop new_id 		
		
	    return scalar boot_0    = mean_risk_b[1]
	    return scalar boot_1    = mean_risk_b[2]
	    return scalar boot_diff = return(boot_1) - return(boot_0)
	    return scalar boot_IPW_95  = IPW_A_95[1]
	    return scalar boot_IPW_99  = IPW_A_99[1]
	    return scalar boot_IPW_max = IPW_A_max[1]
	    
		restore
end		
simulate risk_a0 = r(boot_0) risk_a1 = r(boot_1)	risk_difference=r(boot_diff) IPW_95=r(boot_IPW_95) IPW_99=r(boot_IPW_99)  IPW_max=r(boot_IPW_max), reps(200) seed(12345)  : bootipw_risk_curves

save "$results\Results bootstrap (PP).dta" , replace

***************************

use  "$results\Results bootstrap (PP).dta" , replace

_pctile risk_a0, p(2.5 97.5)
return list

_pctile risk_a1, p(2.5 97.5)
return list

_pctile risk_difference, p(2.5 97.5)
return list

