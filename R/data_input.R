# -------------------------------------------------------------------------------------------

#' Current 'AMECO' data vintage
#'
#' @description Fetches the 'AMECO' data for the EC output gap estimation from the current
#'   vintage.
#'
#' @param country The country name. If left unspecified, data for all countries will be
#'   returned.
#' @param cubs A logical indicating whether the CUBS indicator should be computed if
#'   possible (see details).
#'
#' @details For the computation of CUBS, the following three seasonally adjusted series are
#'   used: the utilization indicators in the service industry, the building and construction
#'   industry, and capacity utilization in manufacturing/industry.
#'
#'  The confidence indicator in the service industry is composed of question 1, 2,
#'   and 3 of the monthly service sector survey  ((Q1 + Q2 + Q3)/3). The underlying survey
#'   questions are as follows:
#'   \itemize{
#'    \item Q1 Business situation development over the past 3 months
#'    \item Q2 Evolution of the demand over the past 3 months
#'    \item Q3 Expectation of the demand over the next 3 months
#'   }
#' The confidence indicator in the building and construction industry is composed of
#'   question 3 and 4 of the monthly building and construction sector survey
#'   ((Q3 and Q4)/2). The underlying survey questions are as follows:
#'   \itemize{
#'    \item Q3 Evolution of your current overall order books
#'    \item Q4 Employment expectations over the next 3 months
#'   }
#' The indicator for capacity utilization in manufacturing/industry is based on question 13
#'   of the quarterly industry sector survey. The underlying survey question is as follows:
#'   \itemize{
#'    \item Q3 Current level of capacity utilization
#'   }
#'
#' @return A list with multiple time series objects for each country. If \code{country} is
#'   specified, a multiple time series object is returned. For each country, the following
#'   series are included:
#'   \item{popw}{Population: 15 to 64 years (unit: 1000 persons, code: NPAN)}
#'   \item{ur}{Unemployment rate, total; Member States: definition EUROSTAT (unit: Percentage of civilian labor force, code: ZUTN)}
#'   \item{etd}{Employment, persons: total economy (National accounts) (unit: 1000 persons, code: NETN)}
#'   \item{et}{Employment, persons: all domestic industries (National accounts) (unit: 1000 persons, code: NETD)}
#'   \item{eet}{Employees, persons: all domestic industries (National accounts) (unit: 1000 persons, code: NWTD)}
#'   \item{vaind}{Gross value added at constant prices: manufacturing industry (unit: bn National currency, code: OVGM)}
#'   \item{vaserv}{Gross value added at constant prices: services (unit: bn National currency, code: OVG5)}
#'   \item{vabuil}{Gross value added at constant prices: building and construction (unit: bn National currency, code: OVG4)}
#'   \item{pconsp}{Price deflator private final consumption expenditure (unit: National currency reference year = 100, code: PCPH)}
#'   \item{cpih}{Harmonised consumer price index (All-items, 2015 = 100, code: ZCPIH)}
#'   \item{cpin}{National consumer price index (All-items, 2015 = 100, code: ZCPIN)}
#'   \item{ngdp}{Gross domestic product at current prices (unit: bn National currency, code: UVGD)}
#'   \item{gdp}{Gross domestic product at constant prices (unit: bn National currency, code: OVGD)}
#'   \item{gdpdefl}{Price deflator gross domestic product (unit: National currency reference year = 100, code: PVGD)}
#'   \item{ahours}{Average annual hours worked per person employed (unit: Hours, code: NLHA)}
#'   \item{l}{Total annual hours worked: total economy (unit: millions, code: NLHT)}
#'   \item{wtotal}{Compensation of employees: total economy (unit: bn National currency, code: UWCD)}
#'   \item{nulc}{Nominal unit labour costs: total economy (Ratio of compensation per employee to real GDP per person employed.) (unit: National currency reference year = 100, code: PLCD)}
#'   \item{k}{Net capital stock at constant prices: total economy (unit: bn National currency, code: OKND)}
#'   \item{serv}{Confidence indicator in the service industry}
#'   \item{buil}{Confidence indicator in the bulding and construction industry}
#'   \item{indu}{Capacity utilization in manufacturing/industry}
#'
#' Additionally, if \code{cubs = TRUE}, the capacity utilization economic sentiment
#'   indicator \code{cubs} will be returned.
#'
#' @source \url{https://economy-finance.ec.europa.eu/economic-research-and-databases/economic-databases/ameco-database/download-annual-data-set-macro-economic-database-ameco_en}
#' @source \url{https://economy-finance.ec.europa.eu/economic-forecast-and-surveys/business-and-consumer-surveys_en}
#' @export
#' @importFrom utils download.file unzip read.delim
fetchAmecoData <- function(country = NULL, cubs = TRUE) {
  
  # temporary fileand folder
  file_path <- tempfile()
  folder <- tempdir()

  # general ameco data -----------------------------------------------

  file_url <- "http://ec.europa.eu/economy_finance/db_indicators/ameco/documents/ameco0.zip"
  # download and unzip
  download.file(url = file_url, destfile = file_path)
  unzip(zipfile = file_path, exdir = folder)
  # file names
  file_names <- paste0(file.path(folder, list.files(path = folder)))
  file_names <- file_names[grepl("AMECO", file_names)]
  
  df <- data.frame()
  for (x in file_names) {
    df_tmp <- read.delim(file = x, sep = ";")
    df <- rbind(df, df_tmp)
  }

  # get list of time series
  tsl <- extract_ameco_data(df = df)

  # indicators for cubs ----------------------------------------------

  # main indicators
  file_url <- "https://ec.europa.eu/economy_finance/db_indicators/surveys/documents/series/nace2_ecfin_2310/main_indicators_sa_nace2.zip"
  file_path <- file.path(folder, paste0("ameco_", Sys.Date(), ".zip"))
  # download and unzip
  download.file(url = file_url, destfile = file_path)
  unzip(zipfile = file_path, exdir = folder)

  # industry total
  file_url <- "https://ec.europa.eu/economy_finance/db_indicators/surveys/documents/series/nace2_ecfin_2310/industry_total_sa_nace2.zip"
  file_path <- file.path(folder, paste0("ameco_", Sys.Date(), ".zip"))
  # download and unzip
  download.file(url = file_url, destfile = file_path)
  unzip(zipfile = file_path, exdir = folder)

  # get list of indicator time series
  tsl_ind <- extract_indicator_data(folder)

  # merge indicators for cubs
  for (x in names(tsl)) {
    namesTmp <- c(colnames(tsl[[x]]), colnames(tsl_ind[[x]]))
    tsl[[x]] <- cbind(tsl[[x]], tsl_ind[[x]])
    colnames(tsl[[x]]) <- namesTmp
  }

  if (cubs) {
    namesCubs <- c("indu", "serv", "buil")
    namesVACubs <- paste0("va", namesCubs)
    for (x in names(tsl)) {
      tryCatch(
        {
          tscubs <- cubs(
            tsCU = tsl[[x]][, namesCubs],
            tsVA = tsl[[x]][, namesVACubs]
          )
          namesTmp <- c(colnames(tsl[[x]]), "cubs")
          tsl[[x]] <- cbind(tsl[[x]], tscubs$cubs)
          colnames(tsl[[x]]) <- namesTmp
        },
        error = function(cond) {  }
      )
    }
  }

  # select country
  if (!is.null(country)) {
    if (country %in% names(tsl)) {
      return(tsl[[country]])
    } else {
      warning("The specified country is not part of the sample. List with all included countries is returned.")
    }
  }

  tsl
}


# -------------------------------------------------------------------------------------------

#' Extracts the relevant 'AMECO' indicator data.
#'
#' @param folder A file path with relevant files.
#'
#' @importFrom zoo as.yearmon as.yearqtr
#' @importFrom openxlsx read.xlsx
#' @keywords internal
extract_indicator_data <- function(folder) {

  # load source and country information
  source <- indicatorList$source
  country <- indicatorList$country

  # load all files
  for (filename in unique(source$file)) {
    sheet <- source$sheet[source$file == filename]
    for (sheetname in sheet) {
      assign(paste0(filename, sheetname), openxlsx::read.xlsx(paste0(folder, "/", filename),
        sheet = sheetname,
        startRow = 1
      ))
    }
  }
  # loop and assign ts object
  for (key in source$key) {
    filename <- source$file[source$key == key]
    sheetname <- source$sheet[source$key == key]
    type <- source$type[source$key == key]
    frequency <- source$freq[source$key == key]
    assign(paste0("key_country", key), sapply(country$abbreviation, function(x) gsub("XXX", x, key)))
    # get start time
    time <- get(paste0(filename, sheetname))[, 1]
    if (frequency == 4) {
      time <- as.yearqtr(time, format = "%Y-Q%q")
      period <- as.numeric(format(time[1], "%q"))
    }
    if (frequency == 12) {
      time <- as.Date(time, origin = "1899-12-30")
      period <- as.numeric(format(time[1], "%m"))
    }
    year <- as.numeric(format(time[1], "%Y"))
    start_date <- c(year, period)
    # ts object
    tmp <- ts(get(paste0(filename, sheetname))[, get(paste0("key_country", key))],
      start = start_date,
      frequency = frequency
    )
    colnames(tmp) <- country$name
    assign(paste0("indicator", type), tmp)
  }
  # reorganize data
  indicator <- list()
  for (c in country$name) {
    tmp <- list()
    for (type in source$type) {
      tmp[[tolower(type)]] <- get(paste0("indicator", type))[, c]
    }
    indicator[[c]] <- tmp
  }
  # temporal aggregation
  tsl_ind <- list()
  for (c in country$name) {
    tmp <- NULL
    mat <- NULL
    for (type in source$type) {
      tmp <- .cubsTa(tsObj = indicator[[c]][[tolower(type)]], conversion = "average", frequency = 1)
      mat <- cbind(mat, tmp)
    }
    colnames(mat) <- tolower(source$type)
    tsl_ind[[c]] <- mat
  }
  tsl_ind
}

# -------------------------------------------------------------------------------------------

#' Extracts the relevant 'AMECO' data
#'
#' @param df A data frame containing all macro-economic 'AMECO' data.
#' @keywords internal
extract_ameco_data <- function(df) {

  # load source and country information
  source <- ameco$source
  country <- ameco$country

  # index and times
  colnames(df) <- gsub("X", "", colnames(df))
  index <- !is.na(as.numeric(colnames(df)))
  years <- as.numeric(colnames(df)[index])
  start <- years[1]
  countrylist <- list()
  # loop and assign ts object
  for (a in country$abbreviation) {
    keynumber <- paste0(a, source$keynumber)
    name <- as.character(country$name[country$abbreviation == a])
    dftmp <- df
    dftmp <- dftmp[dftmp$CODE %in% keynumber, ]
    if (dim(dftmp)[1] > 0) {
      tstmp <- ts(apply(dftmp[, index], 1, function(x) as.numeric(as.character(x))), start = start)
      index_colnames <- sapply(keynumber, match, as.character(dftmp$CODE))
      colnames(tstmp)[index_colnames[!is.na(index_colnames)]] <- as.character(source$finalname[!is.na(index_colnames)])
      countrylist[[name]] <- tstmp
    }
  }

  # delete insufficient data
  countrylist[lapply(countrylist, function(x) dim(x)[2]) < 17] <- NULL
  tsl <- countrylist
  tsl
}


# -------------------------------------------------------------------------------------------

#' Data for estimation
#'
#' @description Computes the necessary input data for the EC output gap estimation on the basis of AMECO
#' data.
#'
#' @param tslAmeco A time series list or a multiple time series object containing AMECO
#'   data.
#' @param alpha A number between \code{0} and \code{1} indicating the labor share. The
#'   default is \code{alpha = 0.65}.
#'
#' @details The list of time series \code{tslAmeco} needs to have the following components:
#' \describe{
#'   \item{popw}{Population: 15 to 64 years (unit: 1000 persons, code: NPAN)}
#'   \item{ur}{Unemployment rate, total; Member States: definition EUROSTAT (unit: Percentage of civilian labor force, code: ZUTN)}
#'   \item{etd}{Employment, persons: total economy (National accounts) (unit: 1000 persons, code: NETN)}
#'   \item{et}{Employment, persons: all domestic industries (National accounts) (unit: 1000 persons, code: NETD)}
#'   \item{eet}{Employees, persons: all domestic industries (National accounts) (unit: 1000 persons, code: NWTD)}
#'   \item{pconsp}{Price deflator private final consumption expenditure (unit: National currency reference year = 100, code: PCPH)}
#'   \item{ngdp}{Gross domestic product at current prices (unit: bn National currency, code: UVGD)}
#'   \item{gdp}{Gross domestic product at constant prices (unit: bn National currency, code: OVGD)}
#'   \item{l}{Total annual hours worked: total economy (unit: millions, code: NLHT)}
#'   \item{wtotal}{Compensation of employees: total economy (unit: bn National currency, code: UWCD)}
#'   \item{nulc}{Nominal unit labour costs: total economy (Ratio of compensation per employee to real GDP per person employed.) (unit: National currency reference year = 100, code: PLCD)}
#'   \item{k}{Net capital stock at constant prices: total economy (unit: bn National currency, code: OKND)}
#' }
#'
#' @return A list of time series containing the same components as the input list \code{tslAmeco}
#'    and the following additional components:
#'   \item{gdpdefl}{Gross domestic product deflator}
#'   \item{tfp}{Total factor productivity}
#'   \item{lfnd}{Labor force non-domestic (unit: 1000 persons)}
#'   \item{parts}{Participation rate}
#'   \item{ahours}{Average hours worked (unit: hours)}
#'   \item{prod}{Labor productivity (unit: real output in millions per person)}
#'   \item{tot}{Terms of trade (pconcp / gdpdefl)}
#'   \item{ws}{Wage share (unit: compensation per unit of nominal output)}
#'   \item{winfl}{Wage inflation}
#'   \item{rulc}{Real unit labor costs}
#'   
#' @export
#' @examples
#' # load data for Germany
#' data("gap")
#' country <- "Germany"
#' tsListRaw <- gap[[country]]
#' tsListInput <- amecoData2input(tslAmeco = tsListRaw)
amecoData2input <- function(tslAmeco, alpha = 0.65) {

  # initialise list
  tsl <- list()
  tslTmp <- list()

  # mts to ts
  tslTmp <- as.list(tslAmeco)

  namesNew <- c("gdp", "ngdp", "k", "et", "ahours", "ur", "wtotal", "gdpdefl", "pconsp", "popw", "l", "etd", "eet", "nulc")
  if (any(!(namesNew %in% names(tslTmp)))) {
    stop(paste0("The variables ", paste0("\"", namesNew[!(namesNew %in% names(tslTmp))], "\"", collapse = ", "), " are missing."))
  }

  tsl[namesNew] <- tslTmp[namesNew]

  ## ----- tfp estimation

  # cubs
  if ("cubs" %in% names(tslTmp)) {
    tsl$cubs <- tslTmp$cubs
  }

  tsl$tfp <- with(tsl, (gdp * 1000) / (l^alpha * (k * 1000)^(1 - alpha)))

  ## ----- labor trend estimation

  # lfnd (labor force non domestic)
  tsl$lfnd <- with(tsl, et - etd)
  
  # parts
  tsl$parts <- with(tsl, etd / (popw * (1 - ur / 100)))
  
  # ahours (unit l: millions, unit et: thousand)
  tsl$ahours <- with(tsl, l / et * 1000)
  
  ## ----- nawru estimation

  # prod
  tsl$prod <- with(tsl, gdp / et)
  
  # gdpdefl
  tsl$gdpdefl <- with(tsl, ngdp / gdp * 100)
  
  # tot
  tsl$tot <- with(tsl, pconsp / gdpdefl)
  
  # ws
  tsl$ws <- with(tsl, wtotal / ngdp)
  
  # infl: wage inflation
  tsl$winfl <- with(tsl, diff(log(wtotal / eet)))
  
  # nulc
  # tsl$nulc <- with(tsl, (wtotal / eet) / (gdp / et))

  # rulc
  tsl$rulc <- with(tsl, nulc / pconsp)
  

  # return list
  return(tsl)
}
