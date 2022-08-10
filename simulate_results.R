library(data.table)
library(splines)
library(lubridate) # for function yday
library(devEMF) # for saving vector images

# -----------------------------------------------------
# get annual data from ONS and estimate long-term trend
# -----------------------------------------------------

# read ONS annual occurence data
# from https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/datasets/deathsrelatedtodrugpoisoningbydateofoccurrence

aod <- read.csv(url('https://raw.githubusercontent.com/danlewer/drd-time-trends/main/ons_occurence_2022.csv'))
setDT(aod)
aod <- aod[order(year)]
aod[, x := .I]

# plot of annual number of drug-related deaths

emf('deathsPerYear.emf', height = 4, width = 5, family = 'Franklin Gothic Book')

par(mar = c(5, 5, 0, 2))
plot(1, type = 'n', xlim = c(1993, 2020), ylim = c(0, 5000), xlab = NA, ylab = NA, axes = F)
axis(1, seq(1995, 2020, 5), pos = 0)
axis(1, c(1993, 2020), pos = 0, labels = F)
axis(2, 0:5 * 1000, formatC(0:5 * 1000, big.mark = ','), pos = 1993, las = 2)
rect(1993, 0, 2020, 5000)
with(aod, {
  points(year, n, pch = 19)
  lines(year, n)
})
title(xlab = 'Year of death', line = 2.5)
title(ylab = 'Number of deaths due to drug poisoning', line = 3.5)

dev.off()

# model with polynomials for time

m1 <- glm(n ~ poly(x, 3), data = aod, family = 'poisson')
m2 <- glm(n ~ poly(x, 4), data = aod, family = 'poisson')
aod[, p1 := predict(m1, newdata = aod, type = 'response')]
aod[, p2 := predict(m2, newdata = aod, type = 'response')]

# model with cubic splines for time

m3 <- glm(n ~ bs(x, knots = seq(5, 25, 5)), data = aod, family = 'poisson')
aod[, p3 := predict(m3, newdata = aod, type = 'response')]

# add fitted lines to plot - cubic splines and polynomials up to 4 are both good

emf('polynomial_or_spline.emf', height = 5, width = 7, family = 'Franklin Gothic Book')

par(mar = c(5, 5, 2, 10), xpd = NA)
with(aod, {
  plot(year, n, type = 'b', ylim = c(2000, 5000), ylab = 'Annual deaths', axes = F, xlab = 'Year')
  lines(year, p1, col = 'red') # polynomials up to 3
  lines(year, p2, col = 'green') # polynomials up to 4
  lines(year, p3, col = 'blue') # cubic splines
})
axis(1, seq(1995, 2020, 5), pos = 2000)
axis(1, c(1993, 2020), pos = 2000, labels = F)
axis(2, 2:5 * 1000, pos = 1993, las = 2)
rect(1993, 2000, 2020, 5000)
ys <- seq(3300, 4000, length.out = 4)
segments(2021, ys, 2024, ys, col = c('black', 'red', 'green', 'blue'))
text(2025, ys, c('Observed data', 'Polynomial 3', 'Polynomial 4', 'Cubic spline'), adj = 0)

dev.off()

# ----------------------
# simulate daily dataset
# ----------------------

dd <- data.table(dt = seq(from = as.Date('1993-01-01', origin = '1970-01-01'), to = as.Date('2018-12-31', origin = '1970-01-01'), by = 'day'))
dd[, yr := year(dt)]
dd[, dayOfYear := yday(dt)]
dd[, x := (yr - 1992) + dayOfYear / 365]
dd[, yr := yr + dayOfYear / 365]
dd[, p1 := predict(m1, newdata = dd, type = 'response')]
dd[, p2 := predict(m2, newdata = dd, type = 'response')]
dd[, p3 := predict(m3, newdata = dd, type = 'response')]

# showing that the results are the same

with(dd[sample(.N, 25)], {
  points(yr, p1, col = 'red')
  points(yr, p2, col = 'green')
  points(yr, p3, col = 'blue')
})

# add time-varying risks
# ----------------------

# weekday

dd[, weekday := weekdays(dt, abbreviate = T)]
wds <- c('Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun')
dd[, weekday := factor(weekday, wds)]
dd[, weekend := weekday == 'Sat' | weekday == 'Sun']

# week number (first week, last week, other)

dd[, dayn := day(dt)]
dd[, mthn := paste0(year(dt), '-', month(dt))]
dd <- dd[, .(daysInMonth = max(dayn)), mthn][dd, on = 'mthn']
dd$weekn <- 'other'
dd$weekn[dd$dayn <= 7] <- 'first'
dd$weekn[dd$dayn >= (dd$daysInMonth - 6)] <- 'last'
dd[, weekn := factor(weekn, c('last', 'other', 'first'))]

# season - assume sin wave in risk, peaking in August

peak <- 213 # 213th day of year, i.e. 1 August 
dd[, season := dayOfYear - peak]
dd[, season := fifelse(season < 0, season + 365, season)]
dd[, season := season / 365 * (pi * 2)]
dd[, season := cos(season)]
dd[, mth := month(dt)]
dd[, mth := factor(mth, 1:12, month.abb)]

# public holidays
# data from Nager.Date: https://date.nager.at/PublicHoliday/Country/GB

publicHolidays <-c('publicholiday.GB.1993.csv', 'publicholiday.GB.1994.csv', 'publicholiday.GB.1995.csv', 'publicholiday.GB.1996.csv', 'publicholiday.GB.1997.csv', 'publicholiday.GB.1998.csv', 'publicholiday.GB.1999.csv', 'publicholiday.GB.2000.csv', 'publicholiday.GB.2001.csv', 'publicholiday.GB.2002.csv', 'publicholiday.GB.2003.csv', 'publicholiday.GB.2004.csv', 'publicholiday.GB.2005.csv', 'publicholiday.GB.2006.csv', 'publicholiday.GB.2007.csv', 'publicholiday.GB.2008.csv', 'publicholiday.GB.2009.csv', 'publicholiday.GB.2010.csv', 'publicholiday.GB.2011.csv', 'publicholiday.GB.2012.csv', 'publicholiday.GB.2013.csv', 'publicholiday.GB.2014.csv', 'publicholiday.GB.2015.csv', 'publicholiday.GB.2016.csv', 'publicholiday.GB.2017.csv', 'publicholiday.GB.2018.csv')
publicHolidays <- paste0('https://raw.githubusercontent.com/danlewer/drd-time-trends/main/public-holidays/', publicHolidays)
publicHolidays <- do.call(rbind, lapply(publicHolidays, read.csv))
publicHolidays <- publicHolidays[, c('Date', 'Name')]
setDT(publicHolidays)
setnames(publicHolidays, c('Date', 'Name'), c('dt', 'publicHoliday'))
publicHolidays[, dt := as.Date(dt)]
dd <- publicHolidays[dd, on = 'dt']
dd[, holiday := !is.na(publicHoliday)]

# create target values risks

weekendRisk <- 1.4
holidayRisk <- 1.3
seasonRisk <- 1.7 # max : min
weekNrisk <- 1.1

dd[, season := (season + 1) / 2 * (seasonRisk - 1) + 1]
dd[, weekendR := fifelse(weekend, weekendRisk, 1)]
dd[, holidayR := fifelse(holiday, holidayRisk, 1)]
dd[, weekNR := fifelse(weekn == 'first', weekNrisk, 1)]
dd[, fRisk := season * weekendR * holidayR * weekNR]
dd[, fRisk2 := exp(log(fRisk) - mean(log(fRisk)))]
dd[, target := (p2 * fRisk2) / 365]

# simulate data
# -------------

set.seed(33)
dd[, sim := rpois(.N, target)]

# -----------------
# fit poisson model
# -----------------

m <- glm(sim ~ poly(x, 4) + weekday + mth + holiday + weekn, data = dd, family = 'poisson')

# -------------------------------------------------------
# report predicted values by weekday, season, and holiday
# -------------------------------------------------------

peak2low <- function(allvars = c('mth', 'weekday', 'holiday', 'weekn'), B = 10, nd) {
  f <- as.formula(paste0('a~poly(x, 4)+', paste0(allvars, collapse = '+')))
  d <- dd[, c('x', allvars), with = F]
  sapply(seq_len(B), function (y) {
    if (y %% 10 == 0) print (y)
    d$a <- rpois(nrow(dd), dd$target)
    m <- glm(f, data = d, family = 'poisson')
    b <- predict(m, newdata = nd, type = 'response')
    max(b) / min(b)
  })
}

critval <- qnorm(0.975)

# example 1: weekday

weekday_predict <- data.table(x = 27, weekday = wds, holiday = F, mth = 'Jan')
weekday_predict_modelled <- predict(m, newdata = weekday_predict, type = 'link', se.fit = T)
weekday_predict[, exp := m$family$linkinv(weekday_predict_modelled$fit)]
weekday_predict[, ul := m$family$linkinv(weekday_predict_modelled$fit + critval * weekday_predict_modelled$se.fit)]
weekday_predict[, ll := m$family$linkinv(weekday_predict_modelled$fit - critval * weekday_predict_modelled$se.fit)]
set.seed(19)
week_peak2low <- peak2low(B = 1000, nd = weekday_predict)
quantile(week_peak2low, probs = c(0.5, 0.025, 0.975))

# example 2: month

month_predict <- data.table(x = 27, mth = month.abb, holiday = F, weekday = 'Mon')
month_predict_modelled <- predict(m, newdata = month_predict, type = 'link', se.fit = T)
month_predict[, exp := m$family$linkinv(month_predict_modelled$fit)]
month_predict[, ul := m$family$linkinv(month_predict_modelled$fit + critval * month_predict_modelled$se.fit)]
month_predict[, ll := m$family$linkinv(month_predict_modelled$fit - critval * month_predict_modelled$se.fit)]
set.seed(17)
month_peak2low <- peak2low(B = 1000, nd = month_predict)
quantile(month_peak2low, probs = c(0.5, 0.025, 0.975))

par(mfrow = c(1, 2))

y <- barplot(weekday_predict$exp, ylim = c(0, 15), main = 'weekday')
arrows(y, weekday_predict$ll, y, weekday_predict$ul, angle = 90, code = 3, length = 0.05)

y <- barplot(month_predict$exp, ylim = c(0, 15), main = 'month')
arrows(y, month_predict$ll, y, month_predict$ul, angle = 90, code = 3, length = 0.05)
