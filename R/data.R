#' Example.A dataset
#'
#'
#' Simulated dataset for demonstrating Full_MEI(), CompareLoadings(), CompareMeans(), and CompareParameters().
#'
#' @format A data frame with 2800 rows and 17 variables:
#' \describe{
#'   \item{ID}{Numeric, case ID number.}
#'   \item{Region}{Character, Grouping variable for Groups: European region ("North West", "Nordic", "Continental", "Southern", "Central East", "North East" and "South East").}
#'   \item{Org Size}{Numeric, control variable - organization size.}
#'   \item{Tenure}{Numeric, control variable - tenure in years.}
#'   \item{R45a}{Numeric, Reverse-coded work-life conflict item 45a: from 1 (never) to 5 (always).}
#'   \item{R45b}{Numeric, Reverse-coded work-life conflict item 45b: from 1 (never) to 5 (always).}
#'   \item{R45c}{Numeric, Reverse-coded work-life conflict item 45c: from 1 (never) to 5 (always).}
#'   \item{R45d}{Numeric, Reverse-coded work-life conflict item 45d: from 1 (never) to 5 (always).}
#'   \item{R45e}{Numeric, Reverse-coded work-life conflict item 45e: from 1 (never) to 5 (always).}
#'   \item{R90a}{Numeric, Reverse-coded engagement item 90a: from 1 (never) to 5 (always).}
#'   \item{R90b}{Numeric, Reverse-coded engagement item 90b: from 1 (never) to 5 (always).}
#'   \item{R90c}{Numeric, Reverse-coded engagement item 90c: from 1 (never) to 5 (always).}
#'   \item{R87a}{Numeric, Reverse-coded psychological wellbeing item 87a: from 1 (at no time) to 6 (all of the time).}
#'   \item{R87b}{Numeric, Reverse-coded psychological wellbeing item 87b: from 1 (at no time) to 6 (all of the time).}
#'   \item{R87c}{Numeric, Reverse-coded psychological wellbeing item 87c: from 1 (at no time) to 6 (all of the time).}
#'   \item{R87d}{Numeric, Reverse-coded psychological wellbeing item 87d: from 1 (at no time) to 6 (all of the time).}
#'   \item{R87e}{Numeric, Reverse-coded psychological wellbeing item 87e: from 1 (at no time) to 6 (all of the time).}
#' }
#' @source Simulated dataset (400 observations for each region) based on the 6th European Working Conditions Survey (EWCS 2015; Eurofound, 2024) -- Examined the life and working conditions of 43,850 respondents in 35 European countries, which were divided into seven regions: North West, Nordic, Continental, Southern, Central East, North East and South East.
"Example.A"


#' Example.B dataset
#'
#'
#' Simulated dataset for demonstrating LGCompareLoadings() and LGCompareMeans() in longitudinal studies.
#'
#' @format A data frame with 242 rows and 16 variables:
#' \describe{
#'   \item{ID}{Numeric, case ID number.}
#'   \item{x1_T1}{Numeric, Satisfaction with life scale Item 1 (SWLS1) at T1 from 1 (lowest satisfaction) to 6 (greatest satisfaction).}
#'   \item{x2_T1}{Numeric, Satisfaction with life scale Item 2 (SWLS1) at T1 from 1 (lowest satisfaction) to 6 (greatest satisfaction).}
#'   \item{x3_T1}{Numeric, Satisfaction with life scale Item 3 (SWLS1) at T1 from 1 (lowest satisfaction) to 6 (greatest satisfaction).}
#'   \item{x4_T1}{Numeric, Satisfaction with life scale Item 4 (SWLS1) at T1 from 1 (lowest satisfaction) to 6 (greatest satisfaction).}
#'   \item{x5_T1}{Numeric, Satisfaction with life scale Item 5 (SWLS1) at T1 from 1 (lowest satisfaction) to 6 (greatest satisfaction).}
#'   \item{x1_T2}{Numeric, Satisfaction with life scale Item 1 (SWLS1) at T2 from 1 (lowest satisfaction) to 6 (greatest satisfaction).}
#'   \item{x2_T2}{Numeric, Satisfaction with life scale Item 2 (SWLS1) at T2 from 1 (lowest satisfaction) to 6 (greatest satisfaction).}
#'   \item{x3_T2}{Numeric, Satisfaction with life scale Item 3 (SWLS1) at T2 from 1 (lowest satisfaction) to 6 (greatest satisfaction).}
#'   \item{x4_T2}{Numeric, Satisfaction with life scale Item 4 (SWLS1) at T2 from 1 (lowest satisfaction) to 6 (greatest satisfaction).}
#'   \item{x5_T2}{Numeric, Satisfaction with life scale Item 5 (SWLS1) at T2 from 1 (lowest satisfaction) to 6 (greatest satisfaction).}
#'   \item{x1_T3}{Numeric, Satisfaction with life scale Item 1 (SWLS1) at T3 from 1 (lowest satisfaction) to 6 (greatest satisfaction).}
#'   \item{x2_T3}{Numeric, Satisfaction with life scale Item 2 (SWLS1) at T3 from 1 (lowest satisfaction) to 6 (greatest satisfaction).}
#'   \item{x3_T3}{Numeric, Satisfaction with life scale Item 3 (SWLS1) at T3 from 1 (lowest satisfaction) to 6 (greatest satisfaction).}
#'   \item{x4_T3}{Numeric, Satisfaction with life scale Item 4 (SWLS1) at T3 from 1 (lowest satisfaction) to 6 (greatest satisfaction).}
#'   \item{x5_T3}{Numeric, Satisfaction with life scale Item 5 (SWLS1) at T3 from 1 (lowest satisfaction) to 6 (greatest satisfaction).}
#' }
#' @source Simulated dataset with 242 observations across three time points based on Study 2 in Wu, C., Chen, L. H., & Tsai, Y. (2009). Longitudinal invariance analysis of the satisfaction with life scale. Personality and Individual Differences, 46, 396-401.
"Example.B"


#' Example.C dataset
#'
#'
#' Simulated dataset for demonstrating LGCompareLoadings() and LGCompareMeans() in congruence studies.
#'
#' @format A data frame with 332 rows and 21 variables:
#' \describe{
#'   \item{ID}{Numeric, case ID number.}
#'   \item{x1_T1}{Numeric, Self-ratings of managerial role item 1  from 1 (not at all effective) to 7 (extremely effective).}
#'   \item{x2_T1}{Numeric, Self-ratings of managerial role item 2  from 1 (not at all effective) to 7 (extremely effective).}
#'   \item{x3_T1}{Numeric, Self-ratings of managerial role item 3  from 1 (not at all effective) to 7 (extremely effective).}
#'   \item{x4_T1}{Numeric, Self-ratings of managerial role item 4  from 1 (not at all effective) to 7 (extremely effective).}
#'   \item{x5_T1}{Numeric, Self-ratings of managerial role item 5  from 1 (not at all effective) to 7 (extremely effective).}
#'   \item{x6_T1}{Numeric, Self-ratings of managerial role item 6  from 1 (not at all effective) to 7 (extremely effective).}
#'   \item{x7_T1}{Numeric, Self-ratings of managerial role item 7  from 1 (not at all effective) to 7 (extremely effective).}
#'   \item{x8_T1}{Numeric, Self-ratings of managerial role item 8  from 1 (not at all effective) to 7 (extremely effective).}
#'   \item{x9_T1}{Numeric, Self-ratings of managerial role item 9  from 1 (not at all effective) to 7 (extremely effective).}
#'   \item{x10_T1}{Numeric, Self-ratings of managerial role item 10  from 1 (not at all effective) to 7 (extremely effective).}
#'   \item{x1_T2}{Numeric, Supervisor-ratings of managerial role item 1  from 1 (not at all effective) to 7 (extremely effective).}
#'   \item{x2_T2}{Numeric, Supervisor-ratings of managerial role item 2  from 1 (not at all effective) to 7 (extremely effective).}
#'   \item{x3_T2}{Numeric, Supervisor-ratings of managerial role item 3  from 1 (not at all effective) to 7 (extremely effective).}
#'   \item{x4_T2}{Numeric, Supervisor-ratings of managerial role item 4  from 1 (not at all effective) to 7 (extremely effective).}
#'   \item{x5_T2}{Numeric, Supervisor-ratings of managerial role item 5  from 1 (not at all effective) to 7 (extremely effective).}
#'   \item{x6_T2}{Numeric, Supervisor-ratings of managerial role item 6  from 1 (not at all effective) to 7 (extremely effective).}
#'   \item{x7_T2}{Numeric, Supervisor-ratings of managerial role item 7  from 1 (not at all effective) to 7 (extremely effective).}
#'   \item{x8_T2}{Numeric, Supervisor-ratings of managerial role item 8  from 1 (not at all effective) to 7 (extremely effective).}
#'   \item{x9_T2}{Numeric, Supervisor-ratings of managerial role item 9  from 1 (not at all effective) to 7 (extremely effective).}
#'   \item{x10_T2}{Numeric, Supervisor-ratings of managerial role item 10  from 1 (not at all effective) to 7 (extremely effective).}
#' }
#' @source Simulated dataset with 332 dyads of mid-level executives and their supervisors based on Ashford, S. J., & Tsui A. S. (1991). Self-regulation for managerial effectiveness: The role of active feedback seeking. Academy of Management Journal, 34, 251-280.
"Example.C"


#' Example.D dataset
#'
#'
#' Simulated dataset for demonstrating MLCompareLoadings() for multi-level confirmatory factor analysis.
#'
#' @format A data frame with 2500 rows and 9 variables:
#' \describe{
#'   \item{ID}{Numeric, Cluster variable}
#'   \item{x1}{Numeric, Negatively worded item 1 of Oldenburg Burnout Inventory from 1 (completely disagree) to 5 (completely agree).}
#'   \item{x2}{Numeric, Negatively worded item 2 of Oldenburg Burnout Inventory from 1 (completely disagree) to 5 (completely agree).}
#'   \item{x3}{Numeric, Negatively worded item 3 of Oldenburg Burnout Inventory from 1 (completely disagree) to 5 (completely agree).}
#'   \item{x4}{Numeric, Negatively worded item 4 of Oldenburg Burnout Inventory from 1 (completely disagree) to 5 (completely agree).}
#'   \item{x5}{Numeric, Negatively worded item 5 of Oldenburg Burnout Inventory from 1 (completely disagree) to 5 (completely agree).}
#'   \item{x6}{Numeric, Negatively worded item 6 of Oldenburg Burnout Inventory from 1 (completely disagree) to 5 (completely agree).}
#'   \item{x7}{Numeric, Negatively worded item 7 of Oldenburg Burnout Inventory from 1 (completely disagree) to 5 (completely agree).}
#'   \item{x8}{Numeric, Negatively worded item 8 of Oldenburg Burnout Inventory from 1 (completely disagree) to 5 (completely agree).}
#' }
#' @source Simulated dataset based on Gruszczynska, E., Basinska, B. A., & Schaufeli, W. B. 2021. Within- and between-person factor structure of the Oldenburg Burnout Inventory: Analysis of a diary study using multilevel confirmatory factor analysis. PLoS ONE 16(5):e0251257. doi: 10.1371/journal.pone.0251257 -- Eight negatively-worded items of the Oldenburg Burnout Inventory (OLBI) to measure burnout of 250 employees for 10 consecutive working days.
"Example.D"
