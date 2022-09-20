library(data.table) # for data reading and manipulation
library(lubridate) # for formatting dates
library(RColorBrewer) # for colour palettes
library(viridisLite) # more colour palettes
library(stringi) # for processing string variables (this script only uses 'stri_trans_totitle')
library(extrafont) # for custom fonts in plots
library(devEMF) # for enhanced metafile graphic device

# ----------------------------------------
# set simulations for peak-to-low estimate
# ----------------------------------------

p2lB <- 10 # 1000 simulations done for results

# ---------
# read data
# ---------

d <- fread("https://raw.githubusercontent.com/danlewer/drd-time-trends/main/ons_dp_opioid.csv")
d$n[is.na(d$n)] <- 0L
d[, o := as.integer(o)]
d[, no := n - o]

#  ------------------
#  add time variables
#  ------------------

# format date

d[, dt := dmy(dt)]
d[, yr := year(dt)]

# days since 1 Jan 1993

d[, day := .I]

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

publicHolidays <-c('publicholiday.GB.1993.csv', 'publicholiday.GB.1994.csv', 'publicholiday.GB.1995.csv', 'publicholiday.GB.1996.csv', 'publicholiday.GB.1997.csv', 'publicholiday.GB.1998.csv', 'publicholiday.GB.1999.csv', 'publicholiday.GB.2000.csv', 'publicholiday.GB.2001.csv', 'publicholiday.GB.2002.csv', 'publicholiday.GB.2003.csv', 'publicholiday.GB.2004.csv', 'publicholiday.GB.2005.csv', 'publicholiday.GB.2006.csv', 'publicholiday.GB.2007.csv', 'publicholiday.GB.2008.csv', 'publicholiday.GB.2009.csv', 'publicholiday.GB.2010.csv', 'publicholiday.GB.2011.csv', 'publicholiday.GB.2012.csv', 'publicholiday.GB.2013.csv', 'publicholiday.GB.2014.csv', 'publicholiday.GB.2015.csv', 'publicholiday.GB.2016.csv', 'publicholiday.GB.2017.csv', 'publicholiday.GB.2018.csv')
publicHolidays <- paste0('https://raw.githubusercontent.com/danlewer/drd-time-trends/main/public-holidays/', publicHolidays)
publicHolidays <- do.call(rbind, lapply(publicHolidays, fread))
publicHolidays <- publicHolidays[, c('Date', 'Name')]
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

# ----------------
# totals of deaths
# ----------------

# total deaths

sum(d$n) # 78583
sum(d$o, na.rm = T) # 32322

# proportion of deaths that are opioid-related, by year

propo <- d[, .(N = sum(n)), yr][d[!is.na(o), .(n = sum(n), o = sum(o), po = sum(o) / sum(n)), yr], on = 'yr']
fwrite(propo, 'by_year.csv')

# ----------
# test model 
# ----------

mTest <- glm(n ~ poly(day, 4) + mth + weekn + weekday + holiday, data = d, family = 'poisson')

# sense-check long-term trend (looks good)

plot(predict(mTest, newdata = data.table(day = 1:nrow(d), mth = 'Jan', weekn = 'last', weekday = 'Mon', holiday = 'no'), type = 'response'), ylab = 'Daily deaths', xlab = 'Days after 1 Jan 1993', type = 'l', ylim = c(0, 15))

#  -------------------
#  create main results
#  -------------------

# function to estimate results

fdl <- function(var = 'mth', 
                outcome = 'n', 
                strat = 1993:2018,
                nd = data.table(weekn = 'first', weekday = 'Mon', holiday = 'no', mth = 'Jan'),
                B = p2lB) {
  dd <- d[yr %in% strat & !is.na(get(outcome))]
  f <- as.formula(paste0(outcome, '~ poly(day, 4) + mth + weekn + weekday + holiday'))
  f2 <- as.formula(paste0(outcome, '~ poly(day, 4) + ', paste(setdiff(c('mth', 'weekn', 'weekday', 'holiday'), var), collapse = '+')))
  model <- glm(f, data = dd, family = 'poisson')
  model2 <- glm(f2, data = dd, family = 'poisson')
  ev <- anova(model, model2, test = 'LRT')$`Pr(>Chi)`[2]
  nd[, which(var == names(nd))] <- NULL
  nd <- cbind(`names<-`(data.table(levels(dd[, get(var)])), var), nd)
  nd$day <- max(dd$day)
  p <- predict(model, newdata = nd, type = 'link', se.fit = T)
  nd$e <- model$family$linkinv(p$fit)
  nd$ul <- model$family$linkinv(p$fit + qnorm(0.975) * p$se.fit)
  nd$ll <- model$family$linkinv(p$fit - qnorm(0.975) * p$se.fit)
  p2la <- NULL # random ratio (max and min can be any levels)
  p2lb <- NULL # fixed ratio (max and min are fixed by point estimate, eg. April vs. September)
  for (i in seq_len(B)) {
    if (i %% 10 == 0) print (i)
    dd$a <- rpois(nrow(dd), dd[, get(outcome)])
    m <- glm(a ~ poly(day, 4) + mth + weekn + weekday + holiday, data = dd, family = 'poisson')
    b <- predict(m, newdata = nd, type = 'response')
    p2la[i] <- max(b) / min(b)
    p2lb[i] <- b[which.max(nd$e)] / b[which.min(nd$e)]
  }
  list(modelled_vals = nd, peak2lowA = p2la, peak2lowB = p2lb, stat_evidence = ev)
}

set.seed(15)
vars <- c('mth', 'weekn', 'weekday', 'holiday')

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


# extract values for table 1

sapply(uD, function (x) x$stat_evidence)
sapply(sD, function (x) x$stat_evidence)

apply(format(round(sapply(uD, function(x) quantile(x$peak2lowB, probs = c(0.5, 0.025, 0.975))), 2), digits = 2, nsmall = 2), 2, function (x) paste0(x[1], ' (', x[2], '-', x[3], ')'))

p2ls <- apply(format(round(sapply(sD, function(x) quantile(x$peak2lowB, probs = c(0.5, 0.025, 0.975))), 2), digits = 2, nsmall = 2), 2, function (x) paste0(x[1], ' (', x[2], '-', x[3], ')'))
fwrite(data.table(p2ls), 'p2ls.csv')

# --------------------
# plot of main results
# --------------------

cols <- data.table(var = c('mth', 'weekn', 'weekday', 'holiday'), col = brewer.pal(4, 'Set2'))
cols <- data.table(var = c('mth', 'weekn', 'weekday', 'holiday'), col = turbo(8)[4:7])


# function for making plot data

mpd <- function (x) {
  p <- lapply(x, function (x) x$modelled_vals)
  p <- mapply(function (x, y) `names<-`(x[, y, with = F], c('lev', 'e', 'ul', 'll')),
              x = p,
              y = lapply(names(p), function (x) c(x, c('e', 'ul', 'll'))),
              SIMPLIFY = F)
  p <- mapply(function (x, y) cbind(var = x, y), 
              x = names(p),
              y = p,
              SIMPLIFY = F)
  p <- rbindlist(p)
  p[, x := var != shift(var)]
  p[, x := fifelse(x == T & !is.na(x), 1.5, 1)]
  p[, x := cumsum(x)]
  cols[p, on = 'var']
}

# make plot data

pD <- mpd(uD)

# make x-axis labels

pD[, lev := stri_trans_totitle(lev)]
pD$lev[pD$lev == 'Ny'] <- 'New Year'
pD$lev[pD$lev == 'No'] <- 'Not a holiday'
pD$lev[pD$lev == 'Other_holiday'] <- 'Other holiday'

png('Fig2.png', height = 6, width = 9, units = 'in', res = 300, family = 'Franklin Gothic Book')

par(xpd = NA, mar = c(5, 4, 0, 7))
plot(1, type = 'n', xlim = c(0, 28), ylim = c(0, 19), axes = F, xlab = NA, ylab = NA)
rect(-0.5, 0, 28, 19)
with(pD, {
  rect(x-1, 0, x, e, col = col)
  arrows(x-0.5, ll, y1 = ul, code = 3, angle = 90, length = 0.05)
})
axis(2, 0:9 * 2, 0:9 * 2, las = 2, pos = -0.5)
axis(1, c(-0.5, 28), labels = F, pos = 0)
text(pD$x-0.5, -0.2, pD$lev, srt = 90, adj = 1)
ys <- seq(10, 19, length.out = 5)
rect(max(pD$x) * 1.05, ys[-length(ys)], max(pD$x) * 1.08, ys[-1], col = rev(cols$col))
text(max(pD$x) * 1.1, ys[-length(ys)] + diff(ys)/2, rev(c('Calendar\nmonth', 'Week of\nmonth', 'Weekday', 'Public\nholiday')), adj = 0)
title(ylab = 'Daily deaths')

dev.off()

# -----------------------------
# plots of stratified variation
# -----------------------------

# general plot function 

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

png('Fig3.png', height = 3, width = 6, units = 'in', res = 300, family = 'Franklin Gothic Book')

par(mfrow = c(1, 3), mar = c(0, 3, 0, 0), xpd = NA, oma = c(4, 2, 0, 0))

sp(x = sD[1:4]$mth, title = '1993-2001', y = 6:11, toff = 0.5)
sp(x = sD[5:8]$mth, title = '2002-2010', y = 5:10, toff = 0.5)
sp(x = sD[9:12]$mth, title = '2011-2018', y = 11:16, toff = 0.5)

mtext('Daily deaths', side = 2, outer = T, cex = 0.7)
mtext('Month of death', side = 1, outer = T, line = 2, cex = 0.7)

dev.off()

emf('FigS1.emf', height = 8, width = 7, family = 'Franklin Gothic Book')

par(mfrow = c(3, 3), mar = c(0, 2, 1.5, 0), xpd = NA, oma = c(4, 3, 0, 0))

sp(x = sD[1:4]$mth, title = '1993-2001', y = 6:10)
title(ylab = 'Daily count of\ndeaths due to drug poisoning')
sp(x = sD[5:8]$mth, title = '2002-2010', y = 5:9)
sp(x = sD[9:12]$mth, title = '2011-2018', y = 11:15)

sp(x = sO[1:4]$mth, title = '1993-2001', y = 3:7)
title(ylab = 'Daily count of\nopioid-related deaths')
sp(x = sO[5:8]$mth, title = '2002-2010', y = 3:7)
sp(x = sO[9:12]$mth, title = '2011-2018', y = 5:9)

sp(x = sN[1:4]$mth, title = '1993-2001', y = 3:7)
title(ylab = 'Daily count of\ndeaths not involving an opioid')
sp(x = sN[5:8]$mth, title = '2002-2010', y = 2:6)
sp(x = sN[9:12]$mth, title = '2011-2018', y = 5:9)

mtext('Month of death', side = 1, outer = T, line = 2, cex = 0.7)

dev.off()

# week of month

emf('FigS2.emf', height = 6.5, width = 5.5, family = 'Franklin Gothic Book')

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

emf('FigS3.emf', height = 8, width = 7, family = 'Franklin Gothic Book')

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

png('Fig4.png', height = 3, width = 6.5, units = 'in', res = 300, family = 'Franklin Gothic Book')

par(mfrow = c(1, 3), mar = c(0, 3, 0, 0), xpd = NA, oma = c(4, 2, 0, 0))

sp(x = uD$holiday, title = 'All drugs', y = 12:20, ln = F, toff = 0.7, xlabs = hols)
sp(x = uO$holiday, title = 'Opioid-related', y = 5:13, ln = F, toff = 0.7, xlabs = hols)
sp(x = uN$holiday, title = 'Not opioid-related', y = 6:14, ln = F, toff = 0.7, xlabs = hols)

mtext('Daily deaths', side = 2, outer = T, cex = 0.7)

dev.off()


emf('FigS4.emf', height = 8, width = 7, family = 'Franklin Gothic Book')

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

d[, mthDay := paste0(mth, '-', day(dt))]
m2 <- glm(n ~ mthDay + poly(day, 4), data = d, family = 'poisson')
nd <- data.table(day = 9496, mthDay = c(paste0('Dec-', 15:31), paste0('Jan-', 1:14)))
pred <- predict(m2, newdata = nd, type = 'link', se.fit = T)
nd$e <- m2$family$linkinv(pred$fit)
nd$ul <- m2$family$linkinv(pred$fit + qnorm(0.975) * pred$se.fit)
nd$ll <- m2$family$linkinv(pred$fit - qnorm(0.975) * pred$se.fit)
nd[, x := .I]
nd[, cols := 'grey90']
nd$cols[nd$mthDay == 'Jan-1'] <- brewer.pal(3, 'Set1')[2]

emf('FigS5.emf', height = 7, width = 10, family = 'Franklin Gothic Book')

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

# Monthly count, smoothed using LOESS

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

png('Fig1.png', height = 5, width = 8, units = 'in', res = 300, family = 'Franklin Gothic Book')

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
