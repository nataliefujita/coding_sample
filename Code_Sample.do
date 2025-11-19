********************************************************************************
* Housekeeping
********************************************************************************

clear all

set more off

global data "~/Documents/DellaVigna/Behavioral_Inequality/Preventive_Health/Data"
global out "~/Documents/DellaVigna/Behavioral_Inequality/Preventive_Health/Results"

********************************************************************************
* NHIS
********************************************************************************

use "$data/Raw/nhis_101925.dta"

keep if inlist(year, 2010, 2013, 2015, 2018)
gen inc_wt = perweight/4
gen wt = sampweight/4
gen obs = 1

gen age_bin = .
	replace age_bin = 1 if inrange(age, 18, 29)
	replace age_bin = 2 if inrange(age, 30, 39)
	replace age_bin = 3 if inrange(age, 40, 49)
	replace age_bin = 4 if inrange(age, 50, 59)
	replace age_bin = 5 if inrange(age, 60, 69)
	replace age_bin = 6 if inrange(age, 70, 99)

gen inc_bin = .
	replace inc_bin = 1 if incfam07on == 11
	replace inc_bin = 2 if incfam07on == 12
	replace inc_bin = 3 if incfam07on == 22
	replace inc_bin = 4 if incfam07on == 23
	replace inc_bin = 5 if incfam07on == 24

gen inc_hh = .
	replace inc_hh = 17499.5 if incfam07on == 11
	replace inc_hh = 42499.5 if incfam07on == 12
	replace inc_hh = 62499.5 if incfam07on == 22
	replace inc_hh = 87499.5 if incfam07on == 23
	replace inc_hh = 150000 if incfam07on == 24
gen rinc_hh = inc_hh * cpi2009 / 1000

* colonoscopy: every 5y for 50+
gen col_pop = inrange(age, 50, 99)
gen col_y = .
	replace col_y = year - colldy if colldy > 0 & colldy < 9996
	replace col_y = 0.5 if collesty == 11
	replace col_y = 1.5 if collesty == 21
	replace col_y = 2.5 if collesty == 31
	replace col_y = 4 if collesty == 41
	replace col_y = 7.5 if collesty == 51
	replace col_y = 15 if collesty == 61
gen col_screen = (col_pop == 1) & (col_y <= 5)

* mammogram: every 2y for 40+
gen mam_pop = (sex == 2) & inrange(age, 40, 99)
gen mam_y = .
	replace mam_y = year - mamldyr if mamldyr > 0 & mamldyr < 9996
	replace mam_y = 0.5 if mamlgyre == 11
	replace mam_y = 1.5 if mamlgyre == 21
	replace mam_y = 2.5 if mamlgyre == 31
	replace mam_y = 4 if mamlgyre == 41
	replace mam_y = 7.5 if mamlgyre == 51
gen mam_screen = (mam_pop == 1) & (mam_y <= 2)

* pap smear: every 3y for 20+
gen pap_pop = (sex == 2) & inrange(age, 20, 99)
gen pap_y = .
	replace pap_y = year - papldyr if papldyr > 0 & papldyr < 9996
	replace pap_y = 0.5 if paplest == 1
	replace pap_y = 1.5 if paplest == 2
	replace pap_y = 2.5 if paplest == 3
	replace pap_y = 4 if paplest == 4
	replace pap_y = 7.5 if paplest == 5
gen pap_screen = (pap_pop == 1) & (pap_y <= 3)

tempfile nhis_data
save `nhis_data'

********************************************************************************
* BRFSS
********************************************************************************

import sasxport5 "$data/Raw/brfss_102825.xpt", clear

rename (_psu _llcpwt _ststr) (psu wt strata)
gen inc_wt = wt
gen obs = 1

gen age_bin = .
	replace age_bin = 1 if inlist(_ageg5yr, 1, 2)
	replace age_bin = 2 if inlist(_ageg5yr, 3, 4)
	replace age_bin = 3 if inlist(_ageg5yr, 5, 6)
	replace age_bin = 4 if inlist(_ageg5yr, 7, 8)
	replace age_bin = 5 if inlist(_ageg5yr, 9, 10)
	replace age_bin = 6 if inlist(_ageg5yr, 11, 12, 13)

gen inc_bin = .
	replace inc_bin = 1 if inrange(income3, 1, 5)
	replace inc_bin = 2 if income3 == 6
	replace inc_bin = 3 if income3 == 7
	replace inc_bin = 4 if income3 == 8
	replace inc_bin = 5 if inrange(income3, 9, 11)

gen inc_hh = .
	replace inc_hh = 4999.5 if income3 == 1
	replace inc_hh = 12499.5 if income3 == 2
	replace inc_hh = 17499.5 if income3 == 3
	replace inc_hh = 22499.5 if income3 == 4
	replace inc_hh = 29.995 if income3 == 5
	replace inc_hh = 42499.5 if income3 == 6
	replace inc_hh = 62499.5 if income3 == 7
	replace inc_hh = 87499.5 if income3 == 8
	replace inc_hh = 124999.5 if income3 == 9
	replace inc_hh = 174999.5 if income3 == 10
	replace inc_hh = 300000 if income3 == 11
gen rinc_hh = inc_hh * 0.684 / 1000

* colonoscopy: every 5y for 50+
gen col_pop = inrange(_ageg5yr, 7, 13)
gen col_screen = (col_pop == 1) & inrange(colntes1, 1, 3)

* mammogram: every 2y for 40+
gen mam_pop = (sexvar == 2) & inrange(_ageg5yr, 5, 13)
gen mam_screen = (mam_pop == 1) & inlist(howlong, 1, 2)

* pap smear: every 3y for 20+
gen pap_pop = (sexvar == 2) & inrange(_ageg5yr, 1, 13)
gen pap_screen = (pap_pop == 1) & inrange(crvclcnc, 1, 3) & (crvclpap == 1)

tempfile brfss_data
save `brfss_data'

********************************************************************************
* Income
********************************************************************************

foreach d in nhis brfss {
	matrix `d'_inc_mat = J(6,8,.)
	use ``d'_data', clear
	qui svyset psu [pw=inc_wt], strata(strata) singleunit(centered)
	qui gen ses = 2
	forvalues x = 1/6 {
		qui svy, subpop(if age_bin == `x'): mean rinc_hh
		qui estat sd
		matrix `d'_inc_mat[`x',1] = r(mean)[1,1]
		matrix `d'_inc_mat[`x',2] = r(sd)[1,1]
		
		qui _pctile rinc_hh if age_bin == `x' [pw=inc_wt], p(10 25 50 75 90)
		matrix `d'_inc_mat[`x',3] = r(r1)
		matrix `d'_inc_mat[`x',4] = r(r2)
		matrix `d'_inc_mat[`x',5] = r(r3)
		matrix `d'_inc_mat[`x',6] = r(r4)
		matrix `d'_inc_mat[`x',7] = r(r5)
		
		qui replace ses = 1 if rinc_hh <= r(r2) & age_bin == `x'
		qui replace ses = 3 if rinc_hh > r(r4) & age_bin == `x'
		
		qui svy, subpop(if age_bin == `x'): total obs
		matrix `d'_inc_mat[`x',8] = r(table)[1,1]
	}
	save ``d'_data', replace
}

frmttable using "$out/nhis_age_inc_tab", ///
statmat(nhis_inc_mat) tex replace fragment ///
sdec(2, 2, 2, 2, 2, 2, 2, 0) ///
ctitles("", "Mean", "SD", "P10", "P25", "P50", "P75", "P90", "N") ///
rtitles("18-29yo" \ "30-39yo" \ "40-49yo" \ "50-59yo" \ "60-69yo" \"70+yo") ///
pretext("\begin{table}[H] \centering \caption{NHIS Income Summary Statistics by Age Group} \vspace*{-5mm}") ///
posttext("\vspace*{-4mm} {\centering \caption*{\begin{footnotesize} \textit{Notes:} This table presents income summary statistics by age group from the 2010, 2013, 2015, and 2018 NHIS. Income refers to total grouped family income and includes all pre-tax income, including hourly wages, salaries, tips, and commissions. Income is imputed as the midpoint for each bin and 1.5 times the topcoded bin. Values are reported in thousands and deflated to 2009 dollars. \end{footnotesize}}} \end{table}")

frmttable using "$out/brfss_age_inc_tab", ///
statmat(brfss_inc_mat) tex replace fragment ///
sdec(2, 2, 2, 2, 2, 2, 2, 0) ///
ctitles("", "Mean", "SD", "P10", "P25", "P50", "P75", "P90", "N") ///
rtitles("18-29yo" \ "30-39yo" \ "40-49yo" \ "50-59yo" \ "60-69yo" \"70+yo") ///
pretext("\begin{table}[H] \centering \caption{BRFSS Income Summary Statistics by Age Group} \vspace*{-5mm}") ///
posttext("\vspace*{-4mm} {\centering \caption*{\begin{footnotesize} \textit{Notes:} This table presents income summary statistics by age group from the 2024 BRFSS. Income refers to annual household income. Income is imputed as the midpoint for each bin and 1.5 times the topcoded bin. Values are reported in thousands and deflated to 2009 dollars. \end{footnotesize}}} \end{table}")

********************************************************************************
* Incidence
********************************************************************************

foreach c in col mam pap {
	foreach d in nhis brfss {
		matrix `c'_`d'_mat = J(3,2,.)
		use ``d'_data', clear
		qui svyset psu [pw=wt], strata(strata) singleunit(centered)
		forvalues x=1/3 {
			qui svy, subpop(if `c'_pop == 1 & ses == `x'): mean `c'_screen
			matrix `c'_`d'_mat[`x',1] = r(table)[1,1]
			
			qui svy, subpop(if `c'_pop == 1 & ses == `x'): total obs
			matrix `c'_`d'_mat[`x', 2] = r(table)[1,1]
		}
		matrix `c'_`d'_mat = J(1,2,.) \ `c'_`d'_mat
	}
	matrix `c'_mat = `c'_nhis_mat, `c'_brfss_mat
}

matrix cancer_mat = col_mat \ mam_mat \ pap_mat

* Cancer Table
frmttable using "$out/cancer_tab", ///
statmat(cancer_mat) tex replace fragment ///
sdec(2, 0, 2, 0) ///
ctitles("", "NHIS", "N", "BRFSS", "N") ///
rtitles("\textbf{Colorectal Cancer}" \ "\hspace{5mm}Low SES" \ "\hspace{5mm}Middle SES" \ "\hspace{5mm}High SES" \ "\textbf{Breast Cancer}" \ "\hspace{5mm}Low SES" \ "\hspace{5mm}Middle SES" \ "\hspace{5mm}High SES" \ "\textbf{Cervical Cancer}" \ "\hspace{5mm}Low SES" \ "\hspace{5mm}Middle SES" \ "\hspace{5mm}High SES") ///
pretext("\begin{table}[H] \centering \caption{Cancer Screenings Incidence} \vspace*{-5mm}") ///
posttext("\vspace*{-4mm} {\centering \caption*{\begin{footnotesize} \textit{Notes:} This table presents participation in cancer screenings by SES group from the 2010, 2013, 2015, and 2018 NHIS and 2024 BRFSS. Low socioeconomic status (SES) is defined by having household income at or below the 25th percentile within each age bracket. High SES is defined by having household income above the 75th percentile within each age bracket. For each survey, we report the share of adults aged 50+ who have received a colonoscopy in the last 5 years, women aged 40+ who have received a mammogram in the last 2 years, and women aged 20+ who have received a Pap smear in the last 3 years in each SES group. N refers to the (weighted) number of adults we consider for a given screening. \end{footnotesize}}} \end{table}")

* Correlation Matrix
foreach d in nhis brfss {
	use ``d'_data', clear
	qui corr col_screen mam_screen pap_screen [aw=wt]
	matrix `d'_corr_mat = r(C)
}

frmttable using "$out/nhis_corr_tab", ///
statmat(nhis_corr_mat) tex replace fragment ///
sdec(2) ///
ctitles("", "Colorectal Cancer", "Breast Cancer", "Cervical Cancer") ///
rtitles("Colorectal Cancer" \ "Breast Cancer" \ "Cervical Cancer") ///
pretext("\begin{table}[H] \centering \caption{Cancer Screenings Correlation} \vspace*{-5mm}") ///
posttext("\vspace*{-4mm} {\centering \caption*{\begin{footnotesize} \textit{Notes:} This table presents the correlation of cancer screenings from the 2010, 2013, 2015, and 2018 NHIS. We consider the share of adults aged 50+ who have received a colonoscopy in the last 5 years, women aged 40+ who have received a mammogram in the last 2 years, and women aged 20+ who have received a Pap smear in the last 3 years. \end{footnotesize}}} \end{table}")
	
frmttable using "$out/brfss_corr_tab", ///
statmat(brfss_corr_mat) tex replace fragment ///
sdec(2) ///
ctitles("", "Colorectal Cancer", "Breast Cancer", "Cervical Cancer") ///
rtitles("Colorectal Cancer" \ "Breast Cancer" \ "Cervical Cancer") ///
pretext("\begin{table}[H] \centering \caption{Cancer Screenings Correlation} \vspace*{-5mm}") ///
posttext("\vspace*{-4mm} {\centering \caption*{\begin{footnotesize} \textit{Notes:} This table presents the correlation of cancer screenings from the 2024 BRFSS. For each survey, we consider the share of adults aged 50+ who have received a colonoscopy in the last 5 years, women aged 40+ who have received a mammogram in the last 2 years, and women aged 20+ who have received a Pap smear in the last 3 years. \end{footnotesize}}} \end{table}")
