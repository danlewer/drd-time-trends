library(data.table)
library(splines)
library(lubridate) # for function yday

# -----------------------------------------------------
# get annual data from ONS and estimate long-term trend
# -----------------------------------------------------

# read ONS annual occurence data
# from https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/datasets/deathsrelatedtodrugpoisoningbydateofoccurrence

aod <- read.csv(url('https://raw.githubusercontent.com/danlewer/drd-time-trends/main/ons_occurence_2022.csv'))
setDT(aod)
aod <- aod[order(year)]
aod[, x := .I]

# model with polynomials for time

m1 <- glm(n ~ poly(x, 3), data = aod, family = 'poisson')
m2 <- glm(n ~ poly(x, 4), data = aod, family = 'poisson')
aod[, p1 := predict(m1, newdata = aod, type = 'response')]
aod[, p2 := predict(m2, newdata = aod, type = 'response')]

# model with cubic splines for time

m3 <- glm(n ~ bs(x, knots = seq(5, 25, 5)), data = aod, family = 'poisson')
aod[, p3 := predict(m3, newdata = aod, type = 'response')]

# add fitted lines to plot - cubic splines and polynomials up to 4 are both good

with(aod, {
  plot(year, n, type = 'b', ylim = c(0, 5000))
  lines(year, p1, col = 'red') # polynomials up to 3
  lines(year, p2, col = 'green') # polynomials up to 4
  lines(year, p3, col = 'blue') # cubic splines
})

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

# add weekly, seasonal, and public holiday risks
# ----------------------------------------------

# weekday

dd[, weekday := weekdays(dt, abbreviate = T)]
wds <- c('Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun')
dd[, weekday := factor(weekday, wds)]
dd[, weekend := weekday == 'Sat' | weekday == 'Sun']

# season - assume sin wave in risk, peaking in August

peak <- 213 # day of year 213, i.e. 1 August 
dd[, season := dayOfYear - peak]
dd[, season := fifelse(season < 0, season + 365, season)]
dd[, season := season / 365 * (pi * 2)]
dd[, season := cos(season)]
dd[, mth := month(dt)]
dd[, mth := factor(mth, 1:12, month.abb)]

# public holidays (just do Christmas for now)

dd[, holiday := dayOfYear %in% c(358:366, 1:2)]

# create target values risks

weekendRisk <- 1.3
holidayRisk <- 1.3
seasonRisk <- 1.3 # max : min

dd[, season := (season + 1) / 2 * (seasonRisk - 1) + 1]
dd[, weekendR := fifelse(weekend, weekendRisk, 1)]
dd[, holidayR := fifelse(holiday, holidayRisk, 1)]
dd[, fRisk := season * weekendR * holidayR]
dd[, fRisk2 := exp(log(fRisk) - mean(log(fRisk)))]
dd[, target := (p2 * fRisk2) / 365]

# simulate data
# -------------

set.seed(33)
dd[, sim := rpois(.N, target)]

# -----------------
# fit poisson model
# -----------------

m <- glm(sim ~ poly(x, 4) + weekday + mth + holiday, data = dd, family = 'poisson')

# -------------------------------------------------------
# report predicted values by weekday, season, and holiday
# -------------------------------------------------------

critval <- qnorm(0.975)

weekday_predict <- data.table(x = 27, weekday = wds, holiday = F, mth = 'Jan')
weekday_predict_modelled <- predict(m, newdata = weekday_predict, type = 'link', se.fit = T)
weekday_predict[, exp := m$family$linkinv(weekday_predict_modelled$fit)]
weekday_predict[, ul := m$family$linkinv(weekday_predict_modelled$fit + critval * weekday_predict_modelled$se.fit)]
weekday_predict[, ll := m$family$linkinv(weekday_predict_modelled$fit - critval * weekday_predict_modelled$se.fit)]

month_predict <- data.table(x = 27, mth = month.abb, holiday = F, weekday = 'Mon')
month_predict_modelled <- predict(m, newdata = month_predict, type = 'link', se.fit = T)
month_predict[, exp := m$family$linkinv(month_predict_modelled$fit)]
month_predict[, ul := m$family$linkinv(month_predict_modelled$fit + critval * month_predict_modelled$se.fit)]
month_predict[, ll := m$family$linkinv(month_predict_modelled$fit - critval * month_predict_modelled$se.fit)]

par(mfrow = c(1, 2))

y <- barplot(weekday_predict$exp, ylim = c(0, 15), main = 'weekday')
arrows(y, weekday_predict$ll, y, weekday_predict$ul, angle = 90, code = 3, length = 0.05)

y <- barplot(month_predict$exp, ylim = c(0, 15), main = 'month')
arrows(y, month_predict$ll, y, month_predict$ul, angle = 90, code = 3, length = 0.05)