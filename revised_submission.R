setwd("~/Q/drd_trends")

library(data.table) # for data reading and manipulation
library(lubridate) # for formatting dates
library(RColorBrewer) # for colour palettes
library(devEMF) # for enhanced metafile graphic device
library(tsModel) # for harmonic terms
library(tseries) # for pacf (partial autocorrelation function)
library(MASS) # for glm.nb
library(lmtest) # for lrtest (comparing poisson with quasipoisson/negative binomial models)

# add transparency to colours
add.alpha <- function(cols, alpha) rgb(t(col2rgb(cols)/255), alpha = alpha)

# ---------
# read data
# ---------

d <- fread("https://raw.githubusercontent.com/danlewer/drd-time-trends/main/ons_dp_opioid.csv")
d$n[is.na(d$n)] <- 0L
d[, o := as.integer(o)]
d[, no := n - o]

# -------------------------------
# compare to poisson distribution
# -------------------------------

mn <- mean(d$n)

rn <- 0:22
dens <- dpois(x = rn, mn)
n <- sapply(rn, function (x) sum(d$n == x))
dens_scaled <- dens * nrow(d)

emf('FigS1.emf', height = 7, width = 7, units = 'in', family = 'Franklin Gothic Book')
par(xpd = NA, mar = c(4, 4, 2, 0))
plot(1, type = 'n', xlim = c(0, 23), ylim = c(0, 1400), xlab = NA, ylab = NA, axes = F)
points(rn, n, pch = 19, col = 'red')
lines(rn, n, col = 'red')
points(rn, dens_scaled, pch = 19)
lines(rn, dens_scaled)
segments(mn, 0, y1 = 1400, lty = 3)
rect(0, 0, 22, 1400)
axis(1, 0:11 * 2, pos = 0)
axis(2, 0:6 * 200, pos = 0, las = 2)
title(xlab = 'Number of deaths per day', line = 2.5)
title(ylab = 'Frequency', line = 2.5)
segments(16, y0 = c(800, 1000), 18, col = c('black', 'red'))
points(c(17, 17), c(800, 1000), pch = 19, col = c('black', 'red'))
text(c(18.2), c(800, 1000), c('Poisson\ndistribution', 'Observed\ndistribution'), adj = 0)
text(mn, 1475, paste0('Mean = ', round(mn, 2)))
dev.off()

#  ------------------
#  add time variables
#  ------------------

# format date

d[, dt := dmy(dt)]
d[, yr := year(dt)]

# days since 1 Jan 1993

d[, day := .I]
d[, day2 := scale(day)]

# weekday

d[, weekday := weekdays(dt, abbreviate = T)]
wds <- c('Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun')
d[, weekday := factor(weekday, wds)]

# week number (first week, last week, other)

d[, dayn := day(dt)]
d[, mthn := paste0(year(dt), '-', month(dt))]
d <- d[, .(daysInMonth = max(dayn)), mthn][d, on = 'mthn']
d$weekn <- 'other'
d$weekn[d$dayn <= 7] <- 'first'
d$weekn[d$dayn >= (d$daysInMonth - 6)] <- 'last'
d[, weekn := factor(weekn, c('first', 'other', 'last'))]

# month

d[, mth := month(dt)]
d[, mth := factor(mth, 1:12, month.abb)]

# public holidays
# data from Nager.Date: https://date.nager.at/PublicHoliday/Country/GB

publicHolidays <- paste0('publicholiday.GB.', 1993:2018, '.csv')
publicHolidays <- paste0('https://raw.githubusercontent.com/danlewer/drd-time-trends/main/public-holidays/', publicHolidays)
publicHolidays <- do.call(rbind, lapply(publicHolidays, fread))
publicHolidays <- publicHolidays[, c('Date', 'Name')]
setnames(publicHolidays, c('Date', 'Name'), c('dt', 'publicHoliday'))
publicHolidays[, dt := as.Date(dt)]
publicHolidays <- unique(publicHolidays)
d <- publicHolidays[d, on = 'dt']
d[, christmas := mth == 'Dec' & day(dt) %in% (24:30)]
d[, ny := (mth == 'Dec' & day(dt) == 31) | (mth == 'Jan' & day(dt) == 1)]
# d[, ny := (mth == 'Jan' & day(dt) == 1)]
d[, holiday := 'no']
d$holiday[!is.na(d$publicHoliday)] <- 'other_holiday'
d$holiday[d$ny] <- 'ny'
d$holiday[d$christmas] <- 'christmas'
d[, holiday := factor(holiday, c('no', 'christmas', 'ny', 'other_holiday'))]

# proportion of the year that has elapsed (to facilitate harmonic terms)

d <- d[, .(nyd = .N), yr][d, on = 'yr']
d[, yd := yday(dt)]
d[, s := yd / nyd] 
d[, s2 := scale(s)]

# periods

d[, period := findInterval(yr, c(1993, 2002, 2011))]
d[, period := factor(period, 1:3, c('1993-2001', '2002-2010', '2011-2018'))]

# ----------------
# totals of deaths
# ----------------

sum(d$n) # 78583
sum(d$o, na.rm = T) # 32322

# --------
# Figure 1
# --------

# Monthly count, smoothed using LOESS

d[, mthDay := paste0(mth, '-', day(dt))]
d[, mthn := cumsum(dayn == 1)]
mt <- d[, .(n = sum(n), jan1 = any(mthDay == 'Jan-1'), yr = max(yr), d = .N), mthn]
mt <- unique(d[, c('mthn', 'mth')])[mt, on = 'mthn']
mm <- glm(n ~ poly(mthn, 4) + mth + offset(log(d)), data = mt, family = 'poisson')
mt[, nd := n / d]
mt[, d := 1]
mt[, mth := 'Jan']
ml <- loess(nd ~ mthn, data = mt, span = 0.05)
mt[, l := predict(ml)]

pred <- predict(mm, newdata = mt, se.fit = T)
mt$e <- mm$family$linkinv(pred$fit)
mt$ul <- mm$family$linkinv(pred$fit + qnorm(0.975) * pred$se.fit)
mt$ll <- mm$family$linkinv(pred$fit - qnorm(0.975) * pred$se.fit)

brk <- 5
brka <- 1
c1 <- brewer.pal(3, 'Set1')[1]

#png('Fig1.png', height = 5, width = 8, units = 'in', res = 300, family = 'Franklin Gothic Book')
cairo_pdf('Fig1.pdf', height = 5, width = 8, family = 'Franklin Gothic Book')

par(mar = c(4, 5, 0, 9), xpd = NA)
plot(1, type = 'n', xlim = range(mt$mthn), ylim = c(0, 15 - brk + brka), axes = F, xlab = NA, ylab = NA)
with(mt[jan1 == T], segments(mthn, 0, y1 = 15 - brk + brka, lty = 3, lwd = 0.5))
with(mt, {
  polygon(c(mthn, rev(mthn)), c(ll, rev(ul)) - brk + brka, col = 'grey90', border = NA)
  points(mthn, nd - brk + brka, pch = 4, cex = 0.5)
  lines(mthn, e - brk + brka)
  lines(mthn, l - brk + brka, col = c1)
})
rect(1, 0, 312, 15 - brk + brka)
with(mt[yr %% 5 == 0 & jan1 == T], {
  axis(1, mthn, labels = F, pos = 0)
  axis(1, mthn, labels = paste0('Jan\n', yr), tick = F, pos = -0.3)
})
axis(2, seq(brk, 15, 1) - brk + brka, labels = seq(brk, 15, 1), pos = 1, las = 2)
axis(2, 0, 0, pos = 1, las = 2)
polygon(x = c(-5, 5, 5, -5), y = c(0.3, 0.5, 0.7, 0.5), col = 'white', border = NA)
segments(x0 = c(-5, -5), y0 = c(0.3, 0.5), x1 = c(5, 5), y1 = c(0.5, 0.7))
title(xlab = 'Date of death', line = 2.5)
title(ylab = 'Daily number of deaths\ndue to drug poisoning')

points(327.5, 10, pch = 4)
rect(320, 7.5, 335, 8.5, col = 'grey90', border = NA)
segments(320, 8, x1 = 335)
segments(320, 6, x1 = 335, col = c1)
text(340, c(6, 8, 10), c('Smoothed\ntrend*', 'Long-term\ntrend (95% CI)', 'Mean daily\ncount per month'), adj = 0)

dev.off()

# -----------
# basic model
# -----------

bm <- glm(n ~ poly(day2, 4, raw = T) + harmonic(s, period = 1, nfreq = 3) + weekn + weekday + holiday, data = d, family = 'poisson')
nd <- data.table(day2 = max(d$day2), weekday = 'Mon', holiday = 'no', weekn = 'other', s = 1:365/365)
nd[, p := predict(bm, newdata = nd, type = 'response')]
max(nd$p) / min(nd$p)

# sense-check long-term trend (looks good)

plot(predict(bm, newdata = data.table(day2 = d$day2, s = d$s, weekn = 'last', weekday = 'Mon', holiday = 'no'), type = 'response'), ylab = 'Daily deaths', xlab = 'Days after 1 Jan 1993', type = 'l', ylim = c(0, 15))

# overdispersion

bm.qp <- glm(n ~ poly(day2, 4, raw = T) + harmonic(s, period = 1, nfreq = 3) + weekn + weekday + holiday, data = d, family = 'quasipoisson')
bm.nb <- glm.nb(n ~ poly(day2, 4, raw = T) + harmonic(s, period = 1, nfreq = 3) + weekn + weekday + holiday, data = d)

lrtest(bm, bm.qp)
lrtest(bm, bm.nb)

# summarise model coefficients

mf <- function (x, digs = 3) {
  y <- cbind(coef(x), confint(x))
  y <- exp(y)
  y <- format(round(y, digs), nsmall = digs, digits = digs)
  y <- paste0(y[,1], '(', y[,2], '-', y[,3], ')')
  y <- gsub('\\(', ' (', gsub(' ', '', y))
  data.frame(var = names(x$coefficients), coef = y)
}

coefs <- list(poisson = `names<-`(mf(bm), c('var', 'poisson')),
              quasi = `names<-`(mf(bm.qp), c('var', 'quasipoisson')),
              negbin = `names<-`(mf(bm.nb), c('var', 'negbin')))
coefs <- Reduce(function(x, y) merge(x, y, by = "var", sort = F), coefs)
fwrite(coefs, 'model_coefficients.csv')

# autocorrelation

emf('FigS2.emf', height = 5, width = 8, units = 'in', family = 'Franklin Gothic Book')

par(mfrow = c(1, 2))
pacf(residuals(bm, type = 'pearson'), lag.max = 365, main = 'Partial autocorrelation function\nof Pearson residuals', xlab = 'Lag (days)')
pacf(residuals(bm, type = 'deviance'), lag.max = 365, main = 'Partial autocorrelation function\nof deviance residuals', xlab = 'Lag (days)')

dev.off()

# ------------------------------------------
# monthly descriptive plots of harmonic terms
# ------------------------------------------

tmp <- d[period == '2011-2018']
m <- glm(n ~ poly(day, 4, raw = T) + harmonic(s, period = 1, nfreq = 3) + weekn + weekday + holiday, data = tmp, family = 'poisson')
p <- predict(m, tmp, se.fit = T)
li <- m$family$linkinv
tmp$fit <- li(p$fit)
tmp$lower <- li(p$fit - p$se.fit * qnorm(0.975))
tmp$upper <- li(p$fit + p$se.fit * qnorm(0.975))
tmp[, dummy := 1]
mt <- tmp[, lapply(.SD, sum), by = mth, .SDcols = c('n', 'fit', 'lower', 'upper', 'dummy')]
mt <- cbind(mth = mt$mth, mt[, lapply(.SD, function (x) x / dummy), .SDcols = c('n', 'fit', 'lower', 'upper')])
mt[, x := .I]

emf('FigS3.emf', height = 5, width = 5, units = 'in', family = 'Franklin Gothic Book')

plot(1, type = 'n', xlim = c(0, 12), ylim = c(8, 11), xlab = NA, ylab = NA, axes = F)
with(mt, {
  rect(x-1, 8, x, n, col = 'grey88')
  lines(x - 0.5, fit, col = 'red')
  polygon(x = c(x - 0.5, rev(x - 0.5)), y = c(lower, rev(upper)), col = add.alpha('red', 0.3), border = NA)
})
axis(2, 8:11, pos = 0, las = 2)
axis(1, 0:12, labels = F, pos = 8)
axis(1, 1:12 - 0.5, month.abb, las = 2, tick = F, pos = 8)
title(ylab = 'Mean daily deaths')

dev.off()

#  -------------------
#  create main results
#  -------------------

# function to estimate results

ndv <- list(weekn = c('first', 'other', 'last'),
            weekday = wds,
            holiday = c('no', 'christmas', 'ny', 'other_holiday'),
            s = 1:365/365)

nd_year <- d[yr == 2018, c('weekn', 'weekday', 'holiday', 's')]
nd_base <- data.table(weekn = 'first', weekday = 'Mon', holiday = 'no', s = 0)

fdl <- function(var = 's',
                outcome = 'n',
                strat = 1993:2018,
                B = 1000) {
  dd <- d[yr %in% strat & !is.na(get(outcome))]
  f <- as.formula(paste0(outcome, '~ poly(day, 4) + harmonic(s, period = 1, nfreq = 3) + weekn + weekday + holiday'))
  f2 <- paste0(outcome, '~ poly(day, 4) + ', paste(setdiff(c('s', 'weekn', 'weekday', 'holiday'), var), collapse = '+'))
  f2 <- gsub(' \\+ s', '+ harmonic(s, period = 1, nfreq = 3)', f2)
  f2 <- as.formula(f2)
  model <- glm(f, data = dd, family = 'poisson')
  model2 <- glm(f2, data = dd, family = 'poisson')
  ev <- anova(model, model2, test = 'LRT')$`Pr(>Chi)`[2]
  
  # year base
  nd_yr <- cbind(nd_year[, var, with = F], nd_base[, -var, with = F])
  nd_yr$day <- max(dd$day)
  nd_yr$e <- predict(model, newdata = nd_yr, type = 'response')
  
  # variable base
  nd_var <- nd_base
  nd_var[, which(var == names(nd_var))] <- NULL
  nd_var <- cbind(`names<-`(data.table(ndv[[var]]), var), nd_var)
  nd_var$day <- max(dd$day)
  
  p <- predict(model, newdata = nd_var, type = 'link', se.fit = T)
  nd_var$e <- model$family$linkinv(p$fit)
  nd_var$ul <- model$family$linkinv(p$fit + qnorm(0.975) * p$se.fit)
  nd_var$ll <- model$family$linkinv(p$fit - qnorm(0.975) * p$se.fit)
  p2la <- NULL # random ratio (max and min can be any levels)
  p2lb <- NULL # fixed ratio (max and min are fixed by point estimate, eg. April vs. September)
  abs_diff <- NULL
  for (i in seq_len(B)) {
    if (i %% 10 == 0) print (i)
    dd$a <- rpois(nrow(dd), dd[, get(outcome)])
    m <- glm(a ~ poly(day, 4) + harmonic(s, period = 1, nfreq = 3) + weekn + weekday + holiday, data = dd, family = 'poisson')
    b <- predict(m, newdata = nd_yr, type = 'response')
    p2la <- cbind(p2la, c(max(b), min(b)))
    p2lb <- cbind(p2lb, c(b[which.max(nd_yr$e)], b[which.min(nd_yr$e)]))
    abs_diff <- cbind(abs_diff, c(sum(b), min(b) * 365))
  }
  list(modelled_vals = nd_var, peak2lowA = p2la, peak2lowB = p2lb, abs_diff = abs_diff, stat_evidence = ev)
}

set.seed(15)
vars <- c('s', 'weekn', 'weekday', 'holiday')

# unstratified results: drug-related deaths

uD <- mapply(fdl,
             var = vars,
             outcome = 'n',
             strat = list(1993:2018),
             SIMPLIFY = F)

# unstratified results: opioid-related deaths

uO <- mapply(fdl,
             var = vars,
             outcome = 'o',
             strat = list(1993:2018),
             SIMPLIFY = F)

# unstratified results: non-opioid-related deaths

uN <- mapply(fdl,
             var = vars,
             outcome = 'no',
             strat = list(1993:2018),
             SIMPLIFY = F)

# stratified drug-related deaths

sD <- mapply(fdl,
             var = rep(vars, 3),
             outcome = 'n',
             strat = rep(c(list(1993:2001), list(2002:2010), list(2011:2018)), each = 4),
             SIMPLIFY = F)

# stratified opioid-related deaths

sO <- mapply(fdl,
             var = rep(vars, 3),
             outcome = 'o',
             strat = rep(c(list(1993:2001), list(2002:2010), list(2011:2018)), each = 4),
             SIMPLIFY = F)

# stratified non-opioid-related deaths

sN <- mapply(fdl,
             var = rep(vars, 3),
             outcome = 'no',
             strat = rep(c(list(1993:2001), list(2002:2010), list(2011:2018)), each = 4),
             SIMPLIFY = F)

# save results

save(uD, uO, uN, sD, sO, sN, file = 'main_results_15jan2023.Rdata')
#load('main_results_15jan2023.Rdata')

# ---------
# main plot
# ---------

x_season <- 1:365/(365/12)
cl <- brewer.pal(4, 'Paired')[3]

#png('main_results.png', height = 5, width = 9.5, units = 'in', res = 300)
cairo_pdf('Fig2.pdf', height = 5, width = 7, family = 'Franklin Gothic Book')

par(xpd = NA, mar = c(7, 4, 0, 0))
plot(1, type = 'n', xlim = c(-1, 31), ylim = c(0, 21), xlab = NA, ylab = NA, axes = F)

with(uD$s$modelled_vals, {
  lines(x_season, e)
  lines(x_season, ll, lty = 3)
  lines(x_season, ul, lty = 3)
})

with(uD$weekn$modelled_vals, {
  rect(13:15 + 0.5, 0, 14:16 + 0.5, e, col = cl)
  arrows(13:15 + 1, ll, y1 = ul, length = 0.05, angle = 90, code = 3)
})

with(uD$weekday$modelled_vals, {
  rect(17:23 + 1, 0, 18:24 + 1, e, col = cl)
  arrows(17:23 + 1.5, ll, y1 = ul, length = 0.05, angle = 90, code = 3)
})

with(uD$holiday$modelled_vals, {
  rect(25:28 + 1.5, 0, 26:29 + 1.5, e, col = cl)
  arrows(25:28 + 2, ll, y1 = ul, length = 0.05, angle = 90, code = 3)
})

axis(2, 0:10 * 2, pos = -1, las = 2)
axis(1, c(0:11, 
          13:16 + 0.5,
          18:25,
          26:30 + 0.5), pos = 0, labels = F)
segments(-1, 0, x1 = 30.5)
text(c(0:11, 14:16, 18:24 + 0.5, 27:30), -1,
     c(paste0('1 ', month.abb), 'First', 'Other', 'Last', wds, 'None', 'Christmas', 'New Year', 'Other'),
     srt = 90, adj = 1)
text(c(6, 15, 21.5, 28.5), -6.5, c('Day of year', 'Week of\nmonth', 'Weekday', 'Holiday'), adj = c(0.5, 1))

title(ylab = 'Mean deaths per day', line = 2)

dev.off()

# --------------------------
# extract values for table 1
# --------------------------

# summary table for main article

extract_pl <- function (x) {
  peak_low_days <- c(which.min(x$s$modelled_vals$e), which.max(x$s$modelled_vals$e))
  peak_low_days <- seq(as.Date('2018-01-01'), as.Date('2018-12-31'), by = 'day')[peak_low_days]
  peak_low <- sapply(x, function (x) x$modelled_vals[,1][c(which.min(x$modelled_vals$e), which.max(x$modelled_vals$e))])
  peak_low[[1]] <- peak_low_days
  sapply(peak_low, function (x) paste0(x[1], ' - ', x[2]))
}

extract_p <- function (x) {
  p_vals <- sapply(x, function (y) y$stat_evidence)
  ifelse(p_vals < 0.001, '<0.001', round(p_vals, 3))
}

extract_p2l <- function (x) {
  p2l <- sapply(x, function (y) quantile(y$peak2lowA[1,] / y$peak2lowA[2,], c(0.5, 0.025, 0.975)))
  p2l <- format(round(p2l, 2), digits = 2, nsmall = 2)
  paste0(p2l[1,], ' (', p2l[2,], '-', p2l[3,], ')')
}

extract_abs <- function (x) {
  adf <- sapply(x, function (y) quantile(y$abs_diff[1,] - y$abs_diff[2,], c(0.5, 0.025, 0.975)))
  adf <- round(adf, 0)
  paste0(adf[1,], ' (', adf[2,], '-', adf[3,], ')')
}

table1 <- data.table(var = vars, cat = extract_pl(uD), evidence = extract_p(uD), p2l = extract_p2l(uD), abs = extract_abs(uD))

fwrite(table1, 'table1.csv')

# stratified table for supplementary

tableS <- lapply(split(sD, rep(1:3, each = 4)), 
                 function (x) data.table(var = vars, 
                                         cat = extract_pl(x), 
                                         evidence = extract_p(x),
                                         p2l = extract_p2l(x), 
                                         abs = extract_abs(x)))

for (i in seq_along(tableS)) {
  tableS[[i]] <- cbind(strata = levels(d$period)[i], tableS[[i]])
}

fwrite(tableS, 'supp_table.csv')

# -----------------------------
# plots of stratified variation
# -----------------------------

# general plot function for season / day of year

ss <- function (x, y = 5:10, title = NA, toff = 0.25) {
  plot(1, type = 'n', xlim = c(1, 365), ylim = range(y), axes = F, xlab = NA, ylab = NA)
  rect(1, min(y), 365, max(y))
  with(x$s$modelled_vals, {
    lines(1:365, e)
    lines(1:365, ll, lty = 3)
    lines(1:365, ul, lty = 3)
  })
  xl <- seq(1, 365, length.out = 13)
  axis(1, xl, pos = min(y), labels = F)
  axis(1, xl[-length(xl)] + diff(xl)/2, month.abb, las = 2, pos = min(y), tick = F)
  axis(2, pos = 1, y, las = 2)
  p <- x$s$stat_evidence
  p <- if(p < 0.001) 'p<0.001' else paste0('p=', format(round(p, 3), digits = 3, nsmall = 3))
  text(183, max(y) - toff, paste0(title, '\n', p))
}

# general plot function for categorical exposures

sp <- function (x = sD[1:4]$mth, title = '1993-2001', y = 6:10, toff = 0.25, ln = T, xlabs = substr(month.abb, 0, 1)) {
  xr <- nrow(x$modelled_vals)
  p <- round(x$stat_evidence, 3)
  p <- if(p == 0) '<0.001' else paste0('=', format(p, nsmall = 3))
  plot(1, type = 'n', xlim = c(0, xr), ylim = range(y), xlab = NA, ylab = NA, axes = F)
  rect(0, min(y), xr, max(y), col = 'white')
  with(x$modelled_vals, {
    points(seq_len(xr) - 0.5, e, pch = 19, cex = 1.2)
    if (ln) lines(seq_len(xr) - 0.5, e)
    arrows(seq_len(xr) - 0.5, ll, y1 = ul, code = 3, angle = 90, length = 0.035)
  })
  axis(2, y, las = 2, pos = 0)
  text(seq_len(xr) - 0.5, min(y) - 0.15, xlabs, adj = c(0.5, 1))
  text(xr/2, max(y) - toff, paste0(title, '\np', p))
}

# season / calendar month

cairo_pdf('Fig3.pdf', height = 3, width = 6, family = 'Franklin Gothic Book')

par(mfrow = c(1, 3), mar = c(0, 3, 0, 0), xpd = NA, oma = c(6, 2, 0, 0))

ss(x = sD[1:4], title = '1993-2001', y = 6:11, toff = 0.5)
ss(x = sD[5:8], title = '2002-2010', y = 5:10, toff = 0.5)
ss(x = sD[9:12], title = '2011-2018', y = 11:16, toff = 0.5)

mtext('Daily deaths', side = 2, outer = T, cex = 0.7)
mtext('Day of death', side = 1, outer = T, line = 3, cex = 0.7)

dev.off()

emf('FigS4.emf', height = 8, width = 7, family = 'Franklin Gothic Book')

par(mfrow = c(3, 3), mar = c(2, 2, 1.5, 0), xpd = NA, oma = c(4, 3, 0, 0))

ss(x = sD[1:4], title = '1993-2001', y = 6:10, toff = 0.35)
title(ylab = 'Daily count of\ndeaths due to drug poisoning')
ss(x = sD[5:8], title = '2002-2010', y = 5:9, toff = 0.35)
ss(x = sD[9:12], title = '2011-2018', y = 11:15, toff = 0.35)

ss(x = sO[1:4], title = '1993-2001', y = 3:7, toff = 0.35)
title(ylab = 'Daily count of\nopioid-related deaths')
ss(x = sO[5:8], title = '2002-2010', y = 3:7, toff = 0.35)
ss(x = sO[9:12], title = '2011-2018', y = 5:9, toff = 0.35)

ss(x = sN[1:4], title = '1993-2001', y = 3:7, toff = 0.35)
title(ylab = 'Daily count of\ndeaths not involving an opioid')
ss(x = sN[5:8], title = '2002-2010', y = 2:6, toff = 0.35)
ss(x = sN[9:12], title = '2011-2018', y = 5:9, toff = 0.35)

mtext('Month of death', side = 1, outer = T, line = 2, cex = 0.7)

dev.off()

# week of month

emf('FigS5.emf', height = 6.5, width = 5.5, family = 'Franklin Gothic Book')

par(mfrow = c(3, 3), mar = c(0, 2, 1.5, 0), xpd = NA, oma = c(4, 3, 0, 0))

sp(x = sD[1:4]$weekn, title = '1993-2001', y = 6:10, ln = F, toff = 0.3, xlabs = c('First', 'Other', 'Last'))
title(ylab = 'Daily count of\ndeaths due to drug poisoning')
sp(x = sD[5:8]$weekn, title = '2002-2010', y = 5:9, ln = F, toff = 0.3, xlabs = c('First', 'Other', 'Last'))
sp(x = sD[9:12]$weekn, title = '2011-2018', y = 11:15, ln = F, toff = 0.3, xlabs = c('First', 'Other', 'Last'))

sp(x = sO[1:4]$weekn, title = '1993-2001', y = 3:7, ln = F, toff = 0.3, xlabs = c('First', 'Other', 'Last'))
title(ylab = 'Daily count of\nopioid-related deaths')
sp(x = sO[5:8]$weekn, title = '2002-2010', y = 3:7, ln = F, toff = 0.3, xlabs = c('First', 'Other', 'Last'))
sp(x = sO[9:12]$weekn, title = '2011-2018', y = 5:9, ln = F, toff = 0.3, xlabs = c('First', 'Other', 'Last'))

sp(x = sN[1:4]$weekn, title = '1993-2001', y = 3:7, ln = F, toff = 0.3, xlabs = c('First', 'Other', 'Last'))
title(ylab = 'Daily count of\ndeaths not involving an opioid')
sp(x = sN[5:8]$weekn, title = '2002-2010', y = 2:6, ln = F, toff = 0.3, xlabs = c('First', 'Other', 'Last'))
sp(x = sN[9:12]$weekn, title = '2011-2018', y = 5:9, ln = F, toff = 0.3, xlabs = c('First', 'Other', 'Last'))

mtext('Week of month', side = 1, outer = T, line = 2, cex = 0.7)

dev.off()

# weekday

emf('FigS6.emf', height = 8, width = 7, family = 'Franklin Gothic Book')

par(mfrow = c(3, 3), mar = c(0, 2, 1.5, 0), xpd = NA, oma = c(4, 3, 0, 0))

sp(x = sD[1:4]$weekday, title = '1993-2001', y = 6:10, ln = T, toff = 0.3, xlabs = wds)
title(ylab = 'Daily count of\ndeaths due to drug poisoning')
sp(x = sD[5:8]$weekday, title = '2002-2010', y = 5:9, ln = T, toff = 0.3, xlabs = wds)
sp(x = sD[9:12]$weekday, title = '2011-2018', y = 11:15, ln = T, toff = 0.3, xlabs = wds)

sp(x = sO[1:4]$weekday, title = '1993-2001', y = 3:7, ln = T, toff = 0.3, xlabs = wds)
title(ylab = 'Daily count of\nopioid-related deaths')
sp(x = sO[5:8]$weekday, title = '2002-2010', y = 3:7, ln = T, toff = 0.3, xlabs = wds)
sp(x = sO[9:12]$weekday, title = '2011-2018', y = 5:9, ln = T, toff = 0.3, xlabs = wds)

sp(x = sN[1:4]$weekday, title = '1993-2001', y = 3:7, ln = T, toff = 0.3, xlabs = wds)
title(ylab = 'Daily count of\ndeaths not involving an opioid')
sp(x = sN[5:8]$weekday, title = '2002-2010', y = 2:6, ln = T, toff = 0.3, xlabs = wds)
sp(x = sN[9:12]$weekday, title = '2011-2018', y = 5:9, ln = T, toff = 0.3, xlabs = wds)

mtext('Weekday', side = 1, outer = T, line = 2, cex = 0.7)

dev.off()

# public holiday

hols <- c('Not\nholiday', 'Christ\n-mas', 'New\nYear', 'Other\nholiday')

#png('Fig4.png', height = 3, width = 6.5, units = 'in', res = 300, family = 'Franklin Gothic Book')
cairo_pdf('Fig4.pdf', height = 3, width = 6.5, family = 'Franklin Gothic Book')

par(mfrow = c(1, 3), mar = c(0, 3, 0, 0), xpd = NA, oma = c(4, 2, 0, 0))

sp(x = uD$holiday, title = 'All drugs', y = 12:20, ln = F, toff = 0.7, xlabs = hols)
sp(x = uO$holiday, title = 'Opioid-related', y = 5:13, ln = F, toff = 0.7, xlabs = hols)
sp(x = uN$holiday, title = 'Not opioid-related', y = 6:14, ln = F, toff = 0.7, xlabs = hols)

mtext('Daily deaths', side = 2, outer = T, cex = 0.7)

dev.off()


emf('FigS7.emf', height = 8, width = 7, family = 'Franklin Gothic Book')

par(mfrow = c(3, 3), mar = c(0, 2, 1.5, 0), xpd = NA, oma = c(4, 3, 0, 0))

sp(x = sD[1:4]$holiday, title = '1993-2001', y = 6:14, ln = F, toff = 0.7, xlabs = hols)
title(ylab = 'Daily count of\ndeaths due to drug poisoning')
sp(x = sD[5:8]$holiday, title = '2002-2010', y = 5:13, ln = F, toff = 0.7, xlabs = hols)
sp(x = sD[9:12]$holiday, title = '2011-2018', y = 11:19, ln = F, toff = 0.7, xlabs = hols)

sp(x = sO[1:4]$holiday, title = '1993-2001', y = 2:10, ln = F, toff = 0.7, xlabs = hols)
title(ylab = 'Daily count of\nopioid-related deaths')
sp(x = sO[5:8]$holiday, title = '2002-2010', y = 3:11, ln = F, toff = 0.7, xlabs = hols)
sp(x = sO[9:12]$holiday, title = '2011-2018', y = 3:11, ln = F, toff = 0.7, xlabs = hols)

sp(x = sN[1:4]$holiday, title = '1993-2001', y = 2:10, ln = F, toff = 0.7, xlabs = hols)
title(ylab = 'Daily count of\ndeaths not involving an opioid')
sp(x = sN[5:8]$holiday, title = '2002-2010', y = 2:10, ln = F, toff = 0.7, xlabs = hols)
sp(x = sN[9:12]$holiday, title = '2011-2018', y = 5:13, ln = F, toff = 0.7, xlabs = hols)

mtext('Public holiday', side = 1, outer = T, line = 2, cex = 0.7)

dev.off()

#  ----------------------------------------
#  exploratory additional analyses / charts
#  ----------------------------------------

# New year effect - showing it is a 1 Jan effect (rather than 31 Dec)

m2 <- glm(n ~ mthDay + poly(day, 4), data = d, family = 'poisson')
nd <- data.table(day = 9496, mthDay = c(paste0('Dec-', 15:31), paste0('Jan-', 1:14)))
pred <- predict(m2, newdata = nd, type = 'link', se.fit = T)
nd$e <- m2$family$linkinv(pred$fit)
nd$ul <- m2$family$linkinv(pred$fit + qnorm(0.975) * pred$se.fit)
nd$ll <- m2$family$linkinv(pred$fit - qnorm(0.975) * pred$se.fit)
nd[, x := .I]
nd[, cols := 'grey90']
nd$cols[nd$mthDay == 'Jan-1'] <- brewer.pal(3, 'Set1')[2]

emf('FigS8.emf', height = 7, width = 10, family = 'Franklin Gothic Book')

par(xpd = NA, mar = c(5, 5, 0, 0))
plot(1, type = 'n', xlim = c(0, 31), ylim = c(10, 21), axes = F, xlab = NA, ylab = NA)
axis(1, 0:31, pos = 10, labels = F)
text(0:30 + 0.5, 9.6, nd$mthDay, srt = 60, adj = 1)
axis(2, 10:21, pos = 0, las = 2)
rect(0, 10, 31, 21, col = 'white')
with(nd, {
  rect(0:30 + 0.07, 10, 1:31 - 0.07, e, col = cols)
  arrows(0:30 + 0.5, ll, y1 = ul, code = 3, angle = 90, length = 0.06)
})
title(ylab = 'Deaths per day')

dev.off()
