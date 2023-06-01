
#' gap data set
#'
#' A dataset containing economic data on various countries from the 'AMECO' 2018 autumn vintage.
#'
#' @format A list object with 53 country time series objects. Each time series object contains 14 variables:
#' \describe{
#'   \item{popw}{Population: 15 to 64 years (unit: 1000 persons, code: NPAN)}
#'   \item{ur}{Unemployment rate, total; Member States: definition EUROSTAT (unit: Percentage of civilian labor force, code: ZUTN)}
#'   \item{etd}{Employment, persons: total economy (National accounts) (unit: 1000 persons, code: NETN)}
#'   \item{et}{Employment, persons: all domestic industries (National accounts) (unit: 1000 persons, code: NETD)}
#'   \item{eet}{Employees, persons: all domestic industries (National accounts) (unit: 1000 persons, code: NWTD)}
#'   \item{vaind}{Gross value added at 2010 prices: manufacturing industry (unit: bn National currency, code: OVGM)}
#'   \item{vaserv}{Gross value added at 2010 prices: services (unit: bn National currency, code: OVG5)}
#'   \item{vabuil}{Gross value added at 2010 prices: building and construction (unit: bn National currency, code: OVG4)}
#'   \item{pconsp}{Price deflator private final consumption expenditure (unit: National currency 2010 = 100, code: PCPH)}
#'   \item{cpih}{Harmonised consumer price index (All-items, 2015 = 100, code: ZCPIH)}
#'   \item{cpin}{National consumer price index (All-items, 2015 = 100, code: ZCPIN)}
#'   \item{ngdp}{Gross domestic product at current prices (unit: bn National currency, code: UVGD)}
#'   \item{gdp}{Gross domestic product at 2010 reference levels (unit: bn National currency, code: OVGD)}
#'   \item{gdpdefl}{Price deflator gross domestic product (unit: National currency 2010 = 100, code: PVGD)}
#'   \item{ahours}{Average annual hours worked per person employed (unit: Hours, code: NLHA)}
#'   \item{l}{Total annual hours worked: total economy (unit: millions, code: NLHT)}
#'   \item{wtotal}{Compensation of employees: total economy (unit: bn National currency, code: UWCD)}
#'   \item{nulc}{Nominal unit labour costs: total economy (Ratio of compensation per employee to real GDP per person employed.) (unit: National currency 2010 = 100, code: PLCD)}
#'   \item{k}{Net capital stock at 2010 prices: total economy (unit: bn National currency, code: OKND)}
#'   \item{serv}{Confidence indicator in the service industry}
#'   \item{buil}{Confidence indicator in the bulding and construction industry}
#'   \item{indu}{Capacity utilization in manufacturing/industry}
#'   }
#' @source \url{https://economy-finance.ec.europa.eu/economic-research-and-databases/economic-databases/ameco-database_en}
"gap"

#' Indicators fo CUBS
#'
#' A dataset containing the service sector confidence indicator, the construction sector confidence indicator and the capacity utilization in manufacturing/industry.
#'
#' A dataset containing the seasonally adjusted utilization indicators in the service industry, the building and construction industry,
#' and capacity utilization in manufacturing/industry for all EU countries and some neighboring countries at different frequencies.
#'
#' The confidence indicator in the service industry is composed of question 1, 2, and 3 of the monthly service sector survey  ((Q1 + Q2 + Q3)/3).
#' The underlying survey questions are as follows:
#'  \itemize{
#'   \item Q1 Business situation development over the past 3 months
#'   \item Q2 Evolution of the demand over the past 3 months
#'   \item Q3 Expectation of the demand over the next 3 months
#' }
#' The confidence indicator in the building and construction industry is composed of question 3 and 4 of the monthly building and construction sector survey  ((Q3 and Q4)/2).
#' The underlying survey questions are as follows:
#'  \itemize{
#'   \item Q3 Evolution of your current overall order books
#'   \item Q4 Employment expectations over the next 3 months
#' }
#' The indicator for capacity utilization in manufacturing/industry is based on question 13 of the quarterly industry sector survey.
#' The underlying survey question is as follows:
#'  \itemize{
#'   \item Q3 Current level of capacity utilization
#' }
#' @format A list with 53 nested country lists with time series objects. Each country list contains 3 time series variables:
#' \describe{
#'   \item{serv}{Confidence indicator in the service industry.}
#'   \item{buil}{Confidence indicator in the bulding and construction industry.}
#'   \item{indu}{Capacity utilization in manufacturing/industry.}
#'   }
#' @source \url{https://economy-finance.ec.europa.eu/economic-forecast-and-surveys/business-and-consumer-surveys_en}
"indicator"
