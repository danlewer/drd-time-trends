library(data.table)
library(lubridate)
library(viridisLite)

setwd("~/drd_time_trends")

d <- fread("ons_drds.csv")
d$n[is.na(d$n)] <- 0L
d[, dt := dmy(dt)]
d[, yr := year(dt)]

# days since 1 Jan 1993

d[, day := .I]

# weekday

d[, weekday := weekdays(dt, abbreviate = T)]
wds <- c('Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun')
d[, weekday := factor(weekday, wds)]
d[, weekend := weekday == 'Sat' | weekday == 'Sun']

# week number (first week, last week, other)

d[, dayn := day(dt)]
d[, mthn := paste0(year(dt), '-', month(dt))]
d <- d[, .(daysInMonth = max(dayn)), mthn][d, on = 'mthn']
d$weekn <- 'other'
d$weekn[d$dayn <= 7] <- 'first'
d$weekn[d$dayn >= (d$daysInMonth - 6)] <- 'last'
d[, weekn := factor(weekn, c('last', 'other', 'first'))]

# month

d[, mth := month(dt)]
d[, mth := factor(mth, 1:12, month.abb)]

# public holidays
# data from Nager.Date: https://date.nager.at/PublicHoliday/Country/GB

publicHolidays <-c('publicholiday.GB.1993.csv', 'publicholiday.GB.1994.csv', 'publicholiday.GB.1995.csv', 'publicholiday.GB.1996.csv', 'publicholiday.GB.1997.csv', 'publicholiday.GB.1998.csv', 'publicholiday.GB.1999.csv', 'publicholiday.GB.2000.csv', 'publicholiday.GB.2001.csv', 'publicholiday.GB.2002.csv', 'publicholiday.GB.2003.csv', 'publicholiday.GB.2004.csv', 'publicholiday.GB.2005.csv', 'publicholiday.GB.2006.csv', 'publicholiday.GB.2007.csv', 'publicholiday.GB.2008.csv', 'publicholiday.GB.2009.csv', 'publicholiday.GB.2010.csv', 'publicholiday.GB.2011.csv', 'publicholiday.GB.2012.csv', 'publicholiday.GB.2013.csv', 'publicholiday.GB.2014.csv', 'publicholiday.GB.2015.csv', 'publicholiday.GB.2016.csv', 'publicholiday.GB.2017.csv', 'publicholiday.GB.2018.csv')
publicHolidays <- paste0('https://raw.githubusercontent.com/danlewer/drd-time-trends/main/public-holidays/', publicHolidays)
publicHolidays <- do.call(rbind, lapply(publicHolidays, read.csv))
publicHolidays <- publicHolidays[, c('Date', 'Name')]
setDT(publicHolidays)
setnames(publicHolidays, c('Date', 'Name'), c('dt', 'publicHoliday'))
publicHolidays[, dt := as.Date(dt)]
d <- publicHolidays[d, on = 'dt']
d[, christmas := mth == 'Dec' & day(dt) %in% (24:30)]
d[, ny := (mth == 'Dec' & day(dt) == 31) | (mth == 'Jan' & day(dt) == 1)]
d[, holiday := 'no']
d$holiday[!is.na(d$publicHoliday)] <- 'other_holiday'
d$holiday[d$ny] <- 'ny'
d$holiday[d$christmas] <- 'christmas'
d[, holiday := factor(holiday, c('no', 'christmas', 'ny', 'other_holiday'))]

# model 

m <- glm(n ~ poly(day, 4) + mth + weekn + weekday + holiday, data = d, family = 'poisson')

anova(m, glm(n ~ poly(day, 4) + weekn + weekday + holiday, data = d, family = 'poisson'), test = 'LRT')
anova(m, glm(n ~ poly(day, 4) + mth + weekday + holiday, data = d, family = 'poisson'), test = 'LRT')
anova(m, glm(n ~ poly(day, 4) + mth + weekn + holiday, data = d, family = 'poisson'), test = 'LRT')
anova(m, glm(n ~ poly(day, 4) + mth + weekn + weekday, data = d, family = 'poisson'), test = 'LRT')

peak2low <- function(allvars = c('mth', 'weekn', 'weekday', 'holiday'), B = 10, nd) {
  f <- as.formula(paste0('a~poly(day, 4)+', paste0(allvars, collapse = '+')))
  dd <- d[, c('day', allvars), with = F]
  sapply(seq_len(B), function (y) {
    if (y %% 10 == 0) print (y)
    dd$a <- rpois(nrow(d), d$n)
    m <- glm(f, data = dd, family = 'poisson')
    b <- predict(m, newdata = nd, type = 'response')
    max(b) / min(b)
  })
}

# function

prd <- function(model = m,
                var = 'weekn', levels = c('first', 'other', 'last'), p2l = 10,
                nd = data.table(day = 9496, weekn = 'first', weekday = 'Mon', holiday = 'no', mth = 'Jan'),
                critval = qnorm(0.975)) {
  nd[, which(var == names(nd))] <- NULL
  nd <- cbind(`names<-`(data.table(levels), var), nd)
  p <- predict(model, newdata = nd, type = 'link', se.fit = T)
  nd$e <- m$family$linkinv(p$fit)
  nd$ul <- m$family$linkinv(p$fit + critval * p$se.fit)
  nd$ll <- m$family$linkinv(p$fit - critval * p$se.fit)
  peak2lowEsts <- if (!is.na(p2l)) {
    peak2low(B = p2l, nd = nd)
  }
  return(list(modelled_values = nd, peak2low = peak2lowEsts))
}

# sense check long term trend

trend_lt <- prd(var = 'day', levels = d$day, p2l = NA)
plot(trend_lt$modelled_values$e)

# weekday

weekday_levs <- c('Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun')
set.seed(19)
trend_weekday <- prd(var = 'weekday', levels = weekday_levs, p2l = 1000)
quantile(trend_weekday$peak2low, probs = c(0.5, 0.025, 0.975))

# month

set.seed(17)
trend_month <- prd(var = 'mth', levels = month.abb, p2l = 1000)
quantile(trend_month$peak2low, probs = c(0.5, 0.025, 0.975))

# week of month

weekn_levs <- c('first', 'other', 'last')
set.seed(99)
trend_weekn <- prd(var = 'weekn', levels = weekn_levs, p2l = 1000)
quantile(trend_weekn$peak2low, probs = c(0.5, 0.025, 0.975))

# holiday

holiday_levs <- c('no', 'other_holiday', 'christmas', 'ny')
set.seed(199)
trend_holiday <- prd(var = 'holiday', levels = holiday_levs, p2l = 1000)
quantile(trend_holiday$peak2low, probs = c(0.5, 0.025, 0.975))

# plots

points <- c(trend_weekday$modelled_values$e,0,
            trend_month$modelled_values$e,0,
            trend_weekn$modelled_values$e,0,
            trend_holiday$modelled_values$e)

ll <- c(trend_weekday$modelled_values$ll,0,
            trend_month$modelled_values$ll,0,
            trend_weekn$modelled_values$ll,0,
            trend_holiday$modelled_values$ll)

ul <- c(trend_weekday$modelled_values$ul,0,
            trend_month$modelled_values$ul,0,
            trend_weekn$modelled_values$ul,0,
            trend_holiday$modelled_values$ul)

cols <- turbo(7)[3:6]
cols2 <- rep(cols, diff(c(0, which(points == 0), length(points))))

sep <- 0.07

png('plot1.png', height = 5, width = 8, units = 'in', res = 300)

par(xpd = NA)
plot(1, type = 'n', xlim = c(0, length(points)), ylim = c(0, 20), axes = F, xlab = NA, ylab = NA)
axis(1, 0:length(points), pos = 0, labels = F)
text(seq_along(points)-0.5, -1, c(weekday_levs, '', month.abb, '', weekn_levs, '', holiday_levs), srt = 60, adj = 1)
axis(2, 0:4 * 5, pos = 0, las = 2)
rect(0, 0, length(points), 20, col = 'grey97')
rect(seq_along(points)-1, 0, seq_along(points), points, col = cols2)
arrows(seq_along(points) - 0.5, ll, y1 = ul, code = 3, angle = 90, length = 0.05)
title(ylab = 'Deaths per day', line = 2)

dev.off()

#  ----------------
#  zoom in on month
#  ================

png('plot2.png', height = 5, width = 5, units = 'in', res = 300)

par(xpd = NA)
plot(1, type = 'n', xlim = c(0, 12), ylim = c(12, 14.5), axes = F, xlab = NA, ylab = NA)
axis(1, 0:12, pos = 0, labels = F)
text(0:11 + 0.5, 11.9, month.abb, srt = 60, adj = 1)
axis(2, seq(12, 14.5, 0.5), pos = 0, las = 2)
rect(0, 12, 12, 14.5, col = 'grey97')
rect(0:11 + sep, 12, 1:12 - sep, trend_month$modelled_values$e, col = cols[2])
arrows(0:11 + 0.5, trend_month$modelled_values$ll, y1 = trend_month$modelled_values$ul, code = 3, angle = 90, length = 0.07)
title(ylab = 'Deaths per day')

dev.off()

#  ---------------------------
#  exploratory New year effect
#  ===========================

d[, mthDay := paste0(mth, '-', day(dt))]
m2 <- glm(n ~ mthDay + poly(day, 4), data = d, family = 'poisson')
nd <- data.table(day = 9496, mthDay = c(paste0('Dec-', 15:31), paste0('Jan-', 1:14)))
pred <- predict(m2, newdata = nd, type = 'link', se.fit = T)
nd$e <- m2$family$linkinv(pred$fit)
critval <- qnorm(0.975)
nd$ul <- m2$family$linkinv(pred$fit + critval * pred$se.fit)
nd$ll <- m2$family$linkinv(pred$fit - critval * pred$se.fit)
nd[, x := .I]
tb <- turbo(12)[3:7]
nd[, cols := tb[5]]
nd$cols[nd$mthDay == 'Jan-1'] <- tb[2]

png('plot3.png', height = 5, width = 7, units = 'in', res = 300)

par(xpd = NA)
plot(1, type = 'n', xlim = c(0, 31), ylim = c(10, 20), axes = F, xlab = NA, ylab = NA)
axis(1, 0:31, pos = 10, labels = F)
text(0:30 + 0.5, 9.6, nd$mthDay, srt = 60, adj = 1)
axis(2, 10:20, pos = 0, las = 2)
rect(0, 10, 31, 20, col = 'grey97')
with(nd, {
  rect(0:30 + sep, 10, 1:31 - sep, e, col = cols)
  arrows(0:30 + 0.5, ll, y1 = ul, code = 3, angle = 90, length = 0.07)
})
title(ylab = 'Deaths per day')

dev.off()

#  ------------------
#  stratified by time
#  ------------------

# peak 2 lows don't work

d[, ys := findInterval(yr, c(1993, 2002, 2011))]
m2 <- lapply(split(d, f = d$ys), function (x) {
  glm(n ~ poly(day, 4) + mth + weekn + weekday + holiday, data = x, family = 'poisson')
})

nd1 <- data.table(day = max(d[ys == 1, day]), weekn = 'first', weekday = 'Mon', holiday = 'no', mth = 'Jan')
nd2 <- data.table(day = max(d[ys == 2, day]), weekn = 'first', weekday = 'Mon', holiday = 'no', mth = 'Jan')
nd3 <- data.table(day = max(d[ys == 3, day]), weekn = 'first', weekday = 'Mon', holiday = 'no', mth = 'Jan')

stratified_results <- lapply(1:3, function (i) {
  mapply(prd,
         model = list(m2[[i]]),
         var = c('weekday', 'mth', 'weekn', 'holiday'),
         levels = list(weekday_levs, month.abb, weekn_levs, holiday_levs),
         nd = list(get(paste0('nd', i))),
         p2l = 1000,
         SIMPLIFY = F)
})

png('plot4.png', height = 9, width = 9, units = 'in', res = 300)

par(mfrow = c(2, 2))

plot(1, type = 'n', ylim = c(0, 15), xlim = c(1, 7), xlab = 'weekday', ylab = 'Deaths per day')
lines(stratified_results[[1]][[1]]$modelled_values$e, type = 'b')
lines(stratified_results[[2]][[1]]$modelled_values$e, type = 'b')
lines(stratified_results[[3]][[1]]$modelled_values$e, type = 'b')

plot(1, type = 'n', ylim = c(0, 15), xlim = c(1, 12), xlab = 'month', ylab = 'Deaths per day')
lines(stratified_results[[1]][[2]]$modelled_values$e, type = 'b')
lines(stratified_results[[2]][[2]]$modelled_values$e, type = 'b')
lines(stratified_results[[3]][[2]]$modelled_values$e, type = 'b')

plot(1, type = 'n', ylim = c(0, 15), xlim = c(1, 3), xlab = 'week-of-month', ylab = 'Deaths per day')
lines(stratified_results[[1]][[3]]$modelled_values$e, type = 'b')
lines(stratified_results[[2]][[3]]$modelled_values$e, type = 'b')
lines(stratified_results[[3]][[3]]$modelled_values$e, type = 'b')

plot(1, type = 'n', ylim = c(0, 20), xlim = c(1, 4), xlab = 'holiday', ylab = 'Deaths per day')
lines(stratified_results[[1]][[4]]$modelled_values$e, type = 'b')
lines(stratified_results[[2]][[4]]$modelled_values$e, type = 'b')
lines(stratified_results[[3]][[4]]$modelled_values$e, type = 'b')

dev.off()

save(stratified_results, file = 'stratified_results_17aug2022.Rdata')

# doesn't work:

quantile(stratified_results[[1]][[2]]$peak2low, probs = c(0.5, 0.975, 0.025))
quantile(stratified_results[[2]][[2]]$peak2low, probs = c(0.5, 0.975, 0.025))
quantile(stratified_results[[3]][[2]]$peak2low, probs = c(0.5, 0.975, 0.025))

# manual peak to low for months

stratified_mth_p2l <- lapply(split(d, f = d$ys), function (x) {
  f <- as.formula('a ~ poly(day, 4) + mth + weekn + weekday + holiday')
  nd <- data.table(day = max(x$day), weekn = 'first', weekday = 'Mon', holiday = 'no', mth = month.abb)
  lapply(seq_len(1000), function (y) {
    if (y %% 50 == 0) print (y)
    x$a <- rpois(nrow(x), x$n)
    m <- glm(f, data = x, family = 'poisson')
    b <- predict(m, newdata = nd, type = 'response')
    list(ratio = max(b) / min(b),
        min_max = month.abb[c(which.min(b), which.max(b))])
  })
})

#sapply(stratified_mth_p2l, quantile, probs = c(0.5, 0.025, 0.975))

s1 <- sapply(stratified_mth_p2l[[1]], function (x) x[[2]])
s2 <- sapply(stratified_mth_p2l[[2]], function (x) x[[2]])
s3 <- sapply(stratified_mth_p2l[[3]], function (x) x[[2]])

fm <- function (x) colSums(sapply(month.abb, function (y) y == x))

mth_profile <- cbind(apply(s1, 1, fm) / 10,
                     apply(s2, 1, fm) / 10,
                     apply(s3, 1, fm) / 10)
colnames(mth_profile) <- c(t(outer(c('1993-2001', '2002-2010', '2011-2018'), c('min', 'max'), paste0)))
