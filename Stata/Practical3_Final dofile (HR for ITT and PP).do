global dataset = "C:\Users\micha\Desktop\Presentations\NASH (2022)\Final\Practical\Final\Final\Final"

use "$dataset\Final Practical TTE (causal survival analysis).dta"  , replace



******************************************************
******************************************************



**************
* QUESTION 1 *
**************

*   check the data

use "$dataset\Final Practical TTE (causal survival analysis).dta"  , replace

desc

* Calculation of the ITT hazard ratios

**************
* QUESTION 2 *
**************

*   Cox regression 

use "$dataset\Final Practical TTE (causal survival analysis).dta"  , replace

* We will use one observation per subject

by id: egen time_Cox=max(actual_time)
by id: egen Y_Cox=max(Y_t1)

keep if time==0

stset time_Cox,f(Y_Cox)

stcox A male age L


***************************************************


**************
* QUESTION 3 *
**************

* Pooled logistic regression

use "$dataset\Final Practical TTE (causal survival analysis).dta"  , replace

drop actual_time
drop if Y_t1==.

* create baseline variables for A and L
gen A_0=A if time==0
gen L_0=L if time==0

by id: egen A_bas=max(A_0)
by id: egen L_bas=max(L_0)


* confounders adjusted in the outcome model
logistic Y_t1 A_bas male age L_bas time c.time#c.time c.time#c.time#c.time, cluster(id)

**************
* QUESTION 4 *
**************

* confounders adjusted through IPW (unstabilised weight)
logistic A_bas male age L_bas if time==0

predict p
gen     pr=p   if A_bas==1
replace pr=1-p if A_bas==0


gen IPW=1/pr

sum pr IPW
logistic Y_t1 A_bas time c.time#c.time c.time#c.time#c.time [pw=IPW], cluster(id)


* confounders adjusted through IPW (stabilised weight)
logistic A_bas male age L_bas if time==0

predict p_den

***

logistic A_bas if time==0

predict p_num

gen     IPW_st=p_num/p_den if A==1
replace IPW_st=(1-p_num)/(1-p_den) if A==0



sum pr IPW IPW_st
logistic Y_t1 A_bas time c.time#c.time c.time#c.time#c.time [pw=IPW_st], cluster(id)



******************************************************
******************************************************

save "$dataset\Practical TTE (before bootstrap).dta" , replace

**************
* QUESTION 6 *
**************

* Calculation of the PP hazard ratios

* First, we need to calculate two variable for deviation from the initial treatment
* One for those who did not adhere from no treatment, while they did not initiate treatment at baseline (A=0 at baseline) --> variable NA_A0
* One for those who did not adhere from treatment, while they initiated treatment at baseline (A=1 at baseline)           --> variable NA_A1

gen NA_A0=0 if A==0 & time==0
gen NA_A1=0 if A==1 & time==0

forvalues i=1(1)4{

replace NA_A0=0 if NA_A0[_n-1]==0 & A==0 & time==`i'
replace NA_A1=0 if NA_A1[_n-1]==0 & A==1 & time==`i'

}

tab time NA_A0, m
tab time NA_A1, m


by id: replace NA_A0=1 if NA_A0==. & NA_A0[_n-1]==0
by id: replace NA_A1=1 if NA_A1==. & NA_A1[_n-1]==0

order id A NA_A0 NA_A1

* we generate a variable for A on the previous time point (i.e. previous year)

by id: gen A_pr=A[_n-1]

order id A A_pr NA_A0 NA_A1 


* First we will calculate the probability of adherence to 
* a) non treamtent and b) treatment
*
* Second, we will then censor the individuals who did not adhere to A=0 and A=1 over time
*
* Third, we will weigh the individuals who actually adhered with the inverse of the probability of adherence to 
* create a pseudo-population where everyone either continuously adheres to treatment or no treatment
*
* Forth, we will estimate the treatment effect from a weighted pooled logistic regression model, 
* strictly in the person-time where individuals are continuously adherent, to get estimates in the pseudo-population.


************************
* Unstabilised weights *
************************


* We calculate the adherence weights for those who did not initiate treatment
* We first calculate the probability of adherence for non treatment initiators

save "$dataset\Practical TTE (analysis adherence).dta" , replace

**** We calculate the Pr(A/A_pr and baseline confounders and time dependent confounders)

* This model runs for time>0 and include years that individuals adhered + one extra year the deviated from the initial strategy

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



order id time A NA_A0 NA_A1  p_NA_A0  p_A_A0  p_A_A1


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

order id A NA_A0 NA_A1  p_NA_A0 p_A_A0 p_A_A1 IPW_A0 IPW_A1 IPW_A



* We now multiply the probabilities at all time points

by id: replace IPW_A=IPW_A*IPW_A[_n-1] if _n~=1

* We now calculate the final (unstabilised) weights 
* by multiplying the adherence weights IPW_A and the weights for randomisation at baseline IPW

gen IPW_final= IPW_A*IPW

sum IPW_final, d


* Few of our weights are rather large (IPW_final>25), however this is <1% of them (67 out of 12872)
* For this reason, we may truncate the few exrteme weights to 25 

* The number 25 has been chosen as a rule of thumb
* Other researchers use as a cutoff for extreme weights from 15 to 50

gen     IPW_final_tr=IPW_final
replace IPW_final_tr=25 if IPW_final>25


save "$dataset\Dataset (before risk curves).dta" , replace


* We can now run the pooled logistic regression using IPW_final_tr

logistic Y_t1 A_bas time c.time#c.time c.time#c.time#c.time [pw=IPW_final_tr], cluster(id)


* Results remain practically the same when we use the non-truncated weights
 
logistic Y_t1 A_bas time c.time#c.time c.time#c.time#c.time [pw=IPW_final], cluster(id)


* or emulate randomisation through outcome regression and weigh the model by IPW_A


logistic Y_t1 A_bas male age L_bas time c.time#c.time c.time#c.time#c.time [pw=IPW_A], cluster(id)

* calculate the truncated weights for IPW_A
gen     IPW_A_tr=IPW_A
replace IPW_A_tr=25 if IPW_A>25

logistic Y_t1 A_bas male age L_bas time c.time#c.time c.time#c.time#c.time [pw=IPW_A_tr], cluster(id)


************************************
************************************
************************************


************************
* Stabilised weights *
************************



* We can also calculate the per-protocol effects using stabilised weights, if our regime is not dynamic


use "$dataset\Practical TTE (analysis adherence).dta" , replace

**** we calculate the Pr(A/A_pr and baseline confounders and time dependent confounders) for those non-treated in the previous year

logistic A male age L_bas L time c.time#c.time c.time#c.time#c.time if A_pr==0 & time>0 & NA_A0~=.

predict p_NA_A0
replace p_NA_A0=. if NA_A0~=0

* probability of non-adherence is set to zero at baseline
replace p_NA_A0=0 if time==0 & A==0

* We now calculate the probability of adherence to no treatment

gen p_A_A0= 1-p_NA_A0


**** we calculate the Pr(A/A_pr and baseline confounders only) , i.e. without L,  for those non-treated in the previous year

logistic A male age L_bas time c.time#c.time c.time#c.time#c.time if A_pr==0 & time>0 & NA_A0~=.

predict p_NA_A0_num
replace p_NA_A0_num=. if NA_A0~=0

* probability of non-adherence is set to zero at baseline
replace p_NA_A0_num=0 if time==0 & A==0

* We now calculate the probability of adherence to no treatment

gen p_A_A0_num= 1-p_NA_A0_num




**** we calculate the Pr(A/A_pr and baseline confounders and time dependent confounders) for those treated in the previous year

logistic A male age L_bas L time c.time#c.time c.time#c.time#c.time  if A_pr==1 & time>0 & NA_A1~=.

predict p_A_A1
replace p_A_A1=. if NA_A1~=0

* probability of adherence is set to one at baseline
replace p_A_A1=1 if time==0 & A==1


**** we calculate the Pr(A/A_pr and baseline confounders only) , i.e. without L,  for those treated in the previous year

logistic A male age L_bas L time c.time#c.time c.time#c.time#c.time  if A_pr==1 & time>0 & NA_A1~=.

predict p_A_A1_num
replace p_A_A1_num=. if NA_A1~=0.

* probability of adherence is set to one at baseline
replace p_A_A1_num=1 if time==0 & A==1




* We drop those who did not adhere (deviated treatment strategy)

drop if NA_A0==1
drop if NA_A1==1
drop if NA_A0==. & NA_A1==.



* We now calculate the stabilised IPW for adherence

* First at each time point
gen IPW_st_A0=p_A_A0_num/p_A_A0
gen IPW_st_A1=p_A_A1_num/p_A_A1


* We now put these weight in one variable

gen     IPW_st_A=IPW_st_A0 if A_bas==0
replace IPW_st_A=IPW_st_A1 if A_bas==1



* We now multiply the probabilities at all time points

by id: replace IPW_st_A=IPW_st_A*IPW_st_A[_n-1] if _n~=1



* We can now run the pooled logistic regression using IPW_st_final

logistic Y_t1 A_bas male age L_bas time c.time#c.time c.time#c.time#c.time [pw=IPW_st_A], cluster(id)




