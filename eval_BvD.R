## Packages
library (boot)
library (car)
library (lme4)
library (blme)
library (DHARMa)
library (pbkrtest)

extract.ci <- function (x) {
    out <- data.frame (numeric (0),
                       numeric (0),
                       numeric (0))
    for (i in 1:length (x [['t0']])) {
        out <- rbind (out, c (x [['t0']] [i],
                              boot.ci (x, index= i,
                                       type= 'perc') [['percent']] [, 4:5]))
    }
    names (out) <- c ('estim', 'lo.ci', 'up.ci')
    out
}


#######################################################################################################
#######################################################################################################
## Characteristics of housing systems
##
## Read data, adjust data types

## Red meat calves
red.df <- read.table ('../Data/inf.red.m.txt', as.is= FALSE, sep= '\t', dec= '.', header= TRUE)
red.df <- red.df [!is.na (red.df [, 'q1.acc.pasture']),]

for (vars in c ('q1.acc.pasture', 'q2.acc.dam', 'q3.suck.dam',
                'q4.transp.6m', 'q5.slaugh.18m')) {
    red.df [red.df [, vars] == 4, vars] <- NA
}
red.df <- red.df [!is.na (red.df [, 'q5.slaugh.18m']),]
red.df <- red.df [!is.na (red.df [, 'q4.transp.6m']),]

for (vars in c ('q1.acc.pasture', 'q2.acc.dam', 'q3.suck.dam',
                'q4.transp.6m', 'q5.slaugh.18m')) {
    red.df [, vars] <- -1 * red.df [, vars] + 4
}

dim (red.df)
names (red.df)
summary (red.df)

# calculate PCA
red.pca <- princomp (red.df [, c ('q1.acc.pasture', 'q2.acc.dam', 'q3.suck.dam',
                                  'q4.transp.6m', 'q5.slaugh.18m')], cor= TRUE)
red.df <- cbind (red.df, predict (red.pca) [, 1:2])
summary (red.pca)
print (loadings (red.pca), cutoff= 0)
plot (red.pca)
biplot (red.pca)
biplot (red.pca, col= c ('white', 'black'))
points (jitter (Comp.2, 20) ~ jitter (Comp.1, 20), red.df,
        col= ifelse (a.origin == 'b', 'red', 'blue'))
plot (jitter (Comp.2, 20) ~ jitter (Comp.1, 20), red.df,
        col= ifelse (a.origin == 'b', 'red', 'blue'))

# Figure
postscript ('Mandel_housing_red.m.eps', width= 5, height= 5,
            horizontal= FALSE, pointsize= 11, paper= 'special')
par (fig= c (0.6, 0.95, 0.55, 0.9), mar= c (0, 0, 0, 0), bty= 'n')
biplot (red.pca, col= c ('white', 'black'),
        xaxt= 'n', yaxt= 'n', xlab= '', ylab= '',
        cex= c (0, 0.5), xlabs= rep ('o', nrow (red.df)),
        ylabs= 1:5)

par (fig= c (0.7, 0.98, 0.15, 0.3), mar= c (0, 4, 0, 0), las= 1, bty= 'n', new= TRUE)
screeplot (red.pca, cex.axis= 0.8, col.lab= 'white', xaxt= 'n', main= '', npcs= 5)
mtext ('Variances', side= 2, line= 2, cex= 0.8, las= 0)

par (fig= c (0, 1, 0, 1), mar= c (3.5, 3.5, 0.2, 0.2), bty= 'o', new= TRUE)
plot (jitter (Comp.2, 20) ~ jitter (Comp.1, 20), red.df,
      col= ifelse (a.origin == 'b', 'red', 'blue'), cex= 1.2,
      xlab= '',
      ylab= '')
mtext ('Extensivity [PC 1]', 1, line= 2.2, las= 0)
mtext ('Early slaughter [PC 2]', 2, line= 2.5, las= 0, cex= 1.1)
dev.off ()


## Replacment calves
replCa.df <- read.table ('../Data/inf.rep.calf.txt', as.is= FALSE, sep= '\t', dec= '.', header= TRUE)
replCa.df <- replCa.df [!is.na (replCa.df [, 'q1.acc.pasture']),]

for (vars in c ('q1.acc.pasture', 'q2.acc.dam', 'q3.suck.dam', 'q4.transp.6m')) {
    replCa.df [replCa.df [, vars] == 4, vars] <- NA
}
replCa.df <- replCa.df [!is.na (replCa.df [, 'q4.transp.6m']),]
replCa.df <- replCa.df [!is.na (replCa.df [, 'q2.acc.dam']),]

for (vars in c ('q1.acc.pasture', 'q2.acc.dam', 'q3.suck.dam', 'q4.transp.6m')) {
    replCa.df [, vars] <- -1 * replCa.df [, vars] + 4
}

dim (replCa.df)
names (replCa.df)
summary (replCa.df)

# calculate PCA
replCa.pca <- princomp (replCa.df [, c ('q1.acc.pasture', 'q2.acc.dam', 'q3.suck.dam', 'q4.transp.6m')],
                        cor= TRUE)
summary (replCa.pca)
replCa.df <- cbind (replCa.df, predict (replCa.pca) [, 1:2])
print (loadings (replCa.pca), cutoff= 0)
plot (replCa.pca)
biplot (replCa.pca)
biplot (replCa.pca, col= c ('white', 'black'))
points (jitter (Comp.2, 20) ~ jitter (Comp.1, 20), replCa.df,
        col= ifelse (a.origin == 'b', 'replCa', 'blue'))
plot (jitter (Comp.2, 20) ~ jitter (Comp.1, 20), replCa.df,
        col= ifelse (a.origin == 'b', 'red', 'blue'))

# Figure
postscript ('Mandel_housing_rep.calf.eps', width= 5, height= 5,
            horizontal= FALSE, pointsize= 11, paper= 'special')
par (fig= c (0.6, 0.95, 0.55, 0.9), mar= c (0, 0, 0, 0), bty= 'n')
biplot (replCa.pca, col= c ('white', 'black'),
        xaxt= 'n', yaxt= 'n', xlab= '', ylab= '',
        cex= c (0, 0.5), xlabs= rep ('o', nrow (replCa.df)),
        ylabs= 1:4)

par (fig= c (0.65, 0.90, 0.33, 0.48), mar= c (0, 4, 0, 0), las= 1, bty= 'n', new= TRUE)
screeplot (replCa.pca, cex.axis= 0.8, col.lab= 'white', xaxt= 'n', main= '', npcs= 4)
mtext ('Variances', side= 2, line= 2, cex= 0.8, las= 0)

par (fig= c (0, 1, 0, 1), mar= c (3.5, 3.5, 0.2, 0.2), bty= 'o', new= TRUE)
plot (jitter (Comp.2, 20) ~ jitter (Comp.1, 20), replCa.df,
      col= ifelse (a.origin == 'b', 'red', 'blue'), cex= 1.2,
      xlab= '',
      ylab= '')
mtext ('Extensivity [PC 1]', 1, line= 2.2, las= 0)
mtext ('Early transport [PC 2]', 2, line= 2.5, las= 0, cex= 1.1)
dev.off ()


## Replacment cows
replCo.df <- read.table ('../Data/inf.rep.cow.txt', as.is= FALSE, sep= '\t', dec= '.', header= TRUE)
replCo.df <- replCo.df [!is.na (replCo.df [, 'q1.acc.pasture']),]

for (vars in c ('q1.acc.pasture', 'q2.indoor', 'q3.slaugh.6y')) {
    replCo.df [replCo.df [, vars] == 4, vars] <- NA
}
replCo.df <- replCo.df [!is.na (replCo.df [, 'q3.slaugh.6y']),]

for (vars in c ('q1.acc.pasture', 'q2.indoor', 'q3.slaugh.6y')) {
    replCo.df [, vars] <- -1 * replCo.df [, vars] + 4
}

dim (replCo.df)
names (replCo.df)
summary (replCo.df)

# calculate PCA
replCo.pca <- princomp (replCo.df [, c ('q1.acc.pasture', 'q2.indoor', 'q3.slaugh.6y')],
                        cor= TRUE)
summary (replCo.pca)
replCo.df <- cbind (replCo.df, predict (replCo.pca) [, 1:2])
print (loadings (replCo.pca), cutoff= 0)
plot (replCo.pca)
biplot (replCo.pca)
biplot (replCo.pca, col= c ('white', 'black'))
points (jitter (Comp.2, 20) ~ jitter (Comp.1, 20), replCo.df,
        col= ifelse (a.origin == 'b', 'red', 'blue'))
plot (jitter (Comp.2, 20) ~ jitter (Comp.1, 20), replCo.df,
        col= ifelse (a.origin == 'b', 'red', 'blue'))

# Figure
postscript ('Mandel_housing_rep.cow.eps', width= 5, height= 5,
            horizontal= FALSE, pointsize= 11, paper= 'special')
par (fig= c (0.64, 0.99, 0.6, 0.95), mar= c (0, 0, 0, 0), bty= 'n')
biplot (replCo.pca, col= c ('white', 'black'),
        xaxt= 'n', yaxt= 'n', xlab= '', ylab= '',
        cex= c (0, 0.5), xlabs= rep ('o', nrow (replCo.df)),
        ylabs= 1:3)

par (fig= c (0.7, 0.95, 0.15, 0.3), mar= c (0, 4, 0, 0), las= 1, bty= 'n', new= TRUE)
screeplot (replCo.pca, cex.axis= 0.8, col.lab= 'white', xaxt= 'n', main= '', npcs= 3)
mtext ('Variances', side= 2, line= 2, cex= 0.8, las= 0)

par (fig= c (0, 1, 0, 1), mar= c (3.5, 3.5, 0.2, 0.2), bty= 'o', new= TRUE)
plot (jitter (Comp.2, 20) ~ jitter (Comp.1, 20), replCo.df,
      col= ifelse (a.origin == 'b', 'red', 'blue'), cex= 1.2,
      xlab= '',
      ylab= '')
mtext ('Access to the outdoors [PC 1]', 1, line= 2.2, las= 0)
mtext ('Early slaughter [PC 2]', 2, line= 2.5, las= 0, cex= 1.1)
dev.off ()


## Veal
veal.df <- read.table ('../Data/inf.veal.txt', as.is= FALSE, sep= '\t', dec= '.', header= TRUE)
veal.df <- veal.df [!is.na (veal.df [, 'q5.slaugh.8m']),]

for (vars in c ('q1.acc.pasture', 'q2.solid.feed', 'q3.bedding',
                'q4.slatted.f', 'q5.slaugh.8m')) {
    veal.df [, vars] <- -1 * veal.df [, vars] + 4
}

dim (veal.df)
names (veal.df)
summary (veal.df)

veal.df <- veal.df [veal.df [, 'a.origin'] == 'd', ]

# calculate PCA
veal.pca <- princomp (veal.df [, c ('q1.acc.pasture', 'q2.solid.feed', 'q3.bedding',
                                    'q4.slatted.f', 'q5.slaugh.8m')],
                        cor= TRUE)
summary (veal.pca)
veal.df <- cbind (veal.df, predict (veal.pca) [, 1:2])
print (loadings (veal.pca), cutoff= 0)
plot (veal.pca)
biplot (veal.pca)
biplot (veal.pca, col= c ('white', 'black'))
points (jitter (Comp.2, 20) ~ jitter (Comp.1, 20), veal.df,
        col= ifelse (a.origin == 'b', 'veal', 'blue'))
plot (jitter (Comp.2, 20) ~ jitter (Comp.1, 20), veal.df,
        col= ifelse (a.origin == 'b', 'red', 'blue'))

# Figure
postscript ('Mandel_housing_veal.eps', width= 5, height= 5,
            horizontal= FALSE, pointsize= 11, paper= 'special')
par (fig= c (0.64, 0.99, 0.65, 0.99), mar= c (0, 0, 0, 0), bty= 'n')
biplot (veal.pca, col= c ('white', 'black'),
        xaxt= 'n', yaxt= 'n', xlab= '', ylab= '',
        cex= c (0, 0.5), xlabs= rep ('o', nrow (veal.df)),
        ylabs= 1:5)

par (fig= c (0.12, 0.4, 0.15, 0.3), mar= c (0, 4, 0, 0), las= 1, bty= 'n', new= TRUE)
screeplot (veal.pca, cex.axis= 0.8, col.lab= 'white', xaxt= 'n', main= '', npcs= 5)
mtext ('Variances', side= 2, line= 2, cex= 0.8, las= 0)

par (fig= c (0, 1, 0, 1), mar= c (3.5, 3.5, 0.2, 0.2), bty= 'o', new= TRUE)
plot (jitter (Comp.2, 20) ~ jitter (Comp.1, 20), veal.df,
      col= ifelse (a.origin == 'b', 'red', 'blue'), cex= 1.2,
      xlab= '',
      ylab= '')
mtext ('Floor comfort [PC 1]', 1, line= 2.2, las= 0)
mtext ('Solid feed, no early slaughter [PC 2]', 2, line= 2.5, las= 0, cex= 1.1)
dev.off ()


#######################################################################################################
#######################################################################################################
## Likelihood-ratings
##
## Read data, adjust data types

meat.df <- read.table ('../Data/aw.rating.txt', as.is= FALSE, sep= '\t', dec= '.', header= TRUE)
meat.df <- meat.df [, 1:23]
meat.df <- meat.df [!is.na (meat.df [, 'q1.food']), ]

dim (meat.df)
names (meat.df)
summary (meat.df)

length (table (meat.df [, 'id']))
table (table (meat.df [, 'id']))
sum (table (table (meat.df [, 'id'])))


## calcualte summary statistics and export as .csv

meat.mean <- aggregate (cbind (q1.food, q2.water, q3.rest.dis, q4.therm.dis,
                               q5.res.mov, q6.injury, q7.disease, q8.pain,
                               q9.norm.behav, q10.other.behav, q11.neg.affect, q12.pos.affect,
                               confidence) ~ origin + categ, meat.df, mean)
meat.sd <- aggregate (cbind (q1.food, q2.water, q3.rest.dis, q4.therm.dis,
                             q5.res.mov, q6.injury, q7.disease, q8.pain,
                             q9.norm.behav, q10.other.behav, q11.neg.affect, q12.pos.affect,
                             confidence) ~ origin + categ, meat.df, sd)
meat.median <- aggregate (cbind (q1.food, q2.water, q3.rest.dis, q4.therm.dis,
                                 q5.res.mov, q6.injury, q7.disease, q8.pain,
                                 q9.norm.behav, q10.other.behav, q11.neg.affect, q12.pos.affect,
                                 confidence) ~ origin + categ, meat.df, median)
meat.min <- aggregate (cbind (q1.food, q2.water, q3.rest.dis, q4.therm.dis,
                              q5.res.mov, q6.injury, q7.disease, q8.pain,
                              q9.norm.behav, q10.other.behav, q11.neg.affect, q12.pos.affect,
                              confidence) ~ origin + categ, meat.df, min)
meat.max <- aggregate (cbind (q1.food, q2.water, q3.rest.dis, q4.therm.dis,
                              q5.res.mov, q6.injury, q7.disease, q8.pain,
                              q9.norm.behav, q10.other.behav, q11.neg.affect, q12.pos.affect,
                              confidence) ~ origin + categ, meat.df, max)
sink ('desc.csv')
for (q in names (meat.mean)[3:15]) {
    for (cat in levels (meat.df [, 'categ'])) {
        for (ori in levels (meat.df [, 'origin'])) {
            cat (formatC (meat.mean [meat.mean [, 'categ'] == cat &
                                     meat.mean [, 'origin'] == ori, q],
                          digits= 1, format= 'f'),
                          ' (',
                 formatC (meat.sd [meat.sd [, 'categ'] == cat &
                                   meat.sd [, 'origin'] == ori, q],
                          digits= 1, format= 'f'),
                 ');', sep= '')
        }
    }
    cat ('\n')
    for (cat in levels (meat.df [, 'categ'])) {
        for (ori in levels (meat.df [, 'origin'])) {
            cat (formatC (meat.median [meat.median [, 'categ'] == cat &
                                       meat.median [, 'origin'] == ori, q],
                          digits= 0, format= 'f'),
                          ' [',
                 formatC (meat.min [meat.min [, 'categ'] == cat &
                                    meat.min [, 'origin'] == ori, q],
                          digits= 0, format= 'f'),
                 '-',
                 formatC (meat.max [meat.max [, 'categ'] == cat &
                                    meat.max [, 'origin'] == ori, q],
                          digits= 0, format= 'f'),
                 '];', sep= '')
        }
    }
    cat ('\n')
}
sink ()



## re-scale likelihood-ratings
#######################################################################################################

# Function:
postscript ('Mandel_trans.eps', width= 5, height= 5, horizontal= FALSE, pointsize= 11, paper= 'special')
par (mar= c (4, 4, 0.5, 0.5), las= 1)
plot (seq (1, 7, 0.5), seq (1/26, 1-1/26, 1/13), type= 'b', pch= 16,
      xlim= c (0.5, 7.5), ylim= c (0, 1),
      xlab= '',
      ylab= 'Re-scaling for logit-transformation')
points (seq (1, 7, 1), seq (1/26, 1-1/26, 2/13), pch= 16, cex= 2)
mtext ('Original Likert scale of the likelihood ratings', side= 1, line= 2.2)

segments (seq (1, 7, 1), rep (0, 7),
          seq (1, 7, 1), seq (1/26, 1-1/26, 2/13), lty= '22')
segments (rep (0.5, 7), seq (1/26, 1-1/26, 2/13),
          seq (1, 7, 1), seq (1/26, 1-1/26, 2/13), lty= '22')
dev.off ()

intToProp <- function (var, scale= c (1, 7), trans= 'logit') {
    cur.seq <- seq (scale [1], scale [2], 0.5)
    vals.range <- length (cur.seq)
    vals.halfdiff <- 1 / (2 * vals.range)
    new.seq <- seq (vals.halfdiff, 1-vals.halfdiff, 2*vals.halfdiff)
    recode.scheme <- data.frame (cur= cur.seq, new= new.seq)
    recode.scheme <- paste (apply (recode.scheme, 1, function (x) paste (x, collapse= '=')),
                            collapse= '; ')
    out <- Recode (var, recode.scheme)

    if (trans == 'logit') out <- logit (out)

    out
}

# Test
intToProp (1:7, trans= 'none')
intToProp (3.5, trans= 'none')
intToProp (1:7)

# on real data on all 12 question complexes
meat.df [, 'q1.food.T'] <- intToProp (meat.df [, 'q1.food'])
meat.df [, 'q2.water.T'] <- intToProp (meat.df [, 'q2.water'])
meat.df [, 'q3.rest.dis.T'] <- intToProp (meat.df [, 'q3.rest.dis'])
meat.df [, 'q4.therm.dis.T'] <- intToProp (meat.df [, 'q4.therm.dis'])
meat.df [, 'q5.res.mov.T'] <- intToProp (meat.df [, 'q5.res.mov'])
meat.df [, 'q6.injury.T'] <- intToProp (meat.df [, 'q6.injury'])
meat.df [, 'q7.disease.T'] <- intToProp (meat.df [, 'q7.disease'])
meat.df [, 'q8.pain.T'] <- intToProp (meat.df [, 'q8.pain'])
meat.df [, 'q9.norm.behav.T'] <- intToProp (meat.df [, 'q9.norm.behav'])
meat.df [, 'q10.other.behav.T'] <- intToProp (meat.df [, 'q10.other.behav'])
meat.df [, 'q11.neg.affect.T'] <- intToProp (meat.df [, 'q11.neg.affect'])
meat.df [, 'q12.pos.affect.T'] <- intToProp (meat.df [, 'q12.pos.affect'])


# prepare variables for PCA
calcResid <- function (var, dat= meat.df) {
    cur.form <- formula (paste (var, ' ~ 1 + (1 | id)', collapse= ''))
    out <- resid (lmer (cur.form, dat))
    plot (qqnorm (out))
    out
}

meat.df [, 'q1.food.R'] <- calcResid ('q1.food.T')
meat.df [, 'q2.water.R'] <- calcResid ('q2.water.T')
meat.df [, 'q3.rest.dis.R'] <- calcResid ('q3.rest.dis.T')
meat.df [, 'q4.therm.dis.R'] <- calcResid ('q4.therm.dis.T')
meat.df [, 'q5.res.mov.R'] <- calcResid ('q5.res.mov.T')
meat.df [, 'q6.injury.R'] <- calcResid ('q6.injury.T')
meat.df [, 'q7.disease.R'] <- calcResid ('q7.disease.T')
meat.df [, 'q8.pain.R'] <- calcResid ('q8.pain.T')
meat.df [, 'q9.norm.behav.R'] <- calcResid ('q9.norm.behav.T')
meat.df [, 'q10.other.behav.R'] <- calcResid ('q10.other.behav.T')
meat.df [, 'q11.neg.affect.R'] <- calcResid ('q11.neg.affect.T')
meat.df [, 'q12.pos.affect.R'] <- calcResid ('q12.pos.affect.T')


## PCA
#######################################################################################################

# calculate PCA
pairs (meat.df [, 36:47])
meat.pca <- princomp (meat.df [, 36:47], cor= TRUE)
summary (meat.pca)
plot (meat.pca)
meat.df <- cbind (meat.df, predict (meat.pca) [, 1:4])

print (loadings (meat.pca), cutoff= 0)
print (loadings (meat.pca), cutoff= 0.4)
biplot (meat.pca, pc.biplot= TRUE)

postscript ('Mandel_pca.eps', width= 5, height= 5, horizontal= FALSE, pointsize= 11, paper= 'special')
par (fig= c (0.5, 0.85, 0.65, 1), mar= c (0, 0, 0, 0), bty= 'n')
biplot (meat.pca, col= c ('white', 'black'),
        xaxt= 'n', yaxt= 'n', xlab= '', ylab= '',
        cex= c (0, 0.5), xlabs= rep ('o', nrow (meat.df)),
        ylabs= 1:12) # paste ('Q', 1:12, sep= ''))

par (fig= c (0.6, 0.95, 0.16, 0.35), mar= c (0, 4, 0, 0), las= 1, bty= 'n', new= TRUE)
screeplot (meat.pca, cex.axis= 0.8, col.lab= 'white', xaxt= 'n', main= '', npcs= 12)
mtext ('Variances', side= 2, line= 1.5, cex= 0.8, las= 0)

par (fig= c (0, 1, 0, 1), mar= c (3.5, 3.5, 0.2, 0.2), bty= 'o', new= TRUE)
plot (Comp.2 ~ Comp.1, meat.df,
      col= ifelse (origin == 'b', 'red', 'blue'), cex= 1.2,
      xlab= '',
      ylab= '')
mtext ('Overall welfare risk [PC 1]', 1, line= 2.2, las= 0)
mtext ('Risk for excessive heat [PC 2]', 2, line= 2.5, las= 0, cex= 1.1)
dev.off ()


## Mixed model
#######################################################################################################

contrasts (meat.df [, 'origin']) <- contr.sum (2)
contrasts (meat.df [, 'categ']) <- contr.sum (4)

meat.df <- cbind (meat.df, model.matrix (~ origin + categ, meat.df) [, -1])
meat.df [, 'oc1'] <- meat.df [, 'origin1'] * meat.df [, 'categ1']
meat.df [, 'oc2'] <- meat.df [, 'origin1'] * meat.df [, 'categ2']
meat.df [, 'oc3'] <- meat.df [, 'origin1'] * meat.df [, 'categ3']
names (meat.df)


# Maximum model
max.lmer <- blmer (Comp.1 ~ origin1 +
                      categ1 + categ2 + categ3 +
                      oc1 + oc2 + oc3 +
                      (1 | id/categ),
                  meat.df, weights= confidence)
summary (max.lmer)

# Residuals
max.simResid <- simulateResiduals (max.lmer, 1000)
plot (max.simResid)
plotResiduals (max.simResid, meat.df [, 'origin'])
plotResiduals (max.simResid, meat.df [, 'categ'])
qqnorm (resid (max.lmer))


# p-values
int.lmer <- blmer (Comp.1 ~ origin1 +
                       categ1 + categ2 + categ3 +
                       (1 | id/categ),
                   meat.df, weights= confidence)
summary (int.lmer)
int.p <- PBmodcomp (max.lmer, int.lmer)
warnings ()
summary (int.p)

BvD.lmer <- blmer (Comp.1 ~ categ1 + categ2 + categ3 +
                       oc1 + oc2 + oc3 +
                       (1 | id/categ),
                   meat.df, weights= confidence)
summary (BvD.lmer)
BvD.p <- PBmodcomp (max.lmer, BvD.lmer)
warnings ()
summary (BvD.p)

AC.lmer <- blmer (Comp.1 ~ origin1 +
                      oc1 + oc2 + oc3 +
                      (1 | id/categ),
                  meat.df, weights= confidence)
summary (AC.lmer)
AC.p <- PBmodcomp (max.lmer, AC.lmer)
warnings ()
summary (AC.p)


# Graphics Component 1
estim.C1.df <- data.frame (origin1= rep (c (1, -1), 4),
                           categ1= rep (c (1, 0, 0, -1), rep (2, 4)),
                           categ2= rep (c (0, 1, 0, -1), rep (2, 4)),
                           categ3= rep (c (0, 0, 1, -1), rep (2, 4)))
estim.C1.df [, 'oc1'] <- estim.C1.df [, 'origin1'] * estim.C1.df [, 'categ1']
estim.C1.df [, 'oc2'] <- estim.C1.df [, 'origin1'] * estim.C1.df [, 'categ2']
estim.C1.df [, 'oc3'] <- estim.C1.df [, 'origin1'] * estim.C1.df [, 'categ3']

estim.C1.Mod <- function (x) predict (x, estim.C1.df, re.form= NA)
estim.C1.raw <- bootMer (max.lmer, estim.C1.Mod, nsim= 1000, .progress= 'win')
estim.C1.val <- extract.ci (estim.C1.raw)

postscript ('Mandel_C1.eps', width= 7, height= 8, horizontal= FALSE, pointsize= 14, paper= 'special')
par (bty= 'n', las= 1, mar= c (4.5, 4, 0.1, 0.2))
boxplot (Comp.1 ~ origin + categ, meat.df,
         ylab= '', xlab= '', xaxt= 'n', yaxt= 'n',
         xlim= c (4/7, 9.5 + 3/7), at= c (1, 2, 3.5, 4.5, 6, 7, 8.5, 9.5),
         range= 0, col= c ('indianred', 'lightblue'), whisklty= 1,
         whiskcol= c ('red3', 'blue3'),
         boxcol= c ('red3', 'blue3'),
         medcol= c ('red3', 'blue3'),
         staplecol= c ('red3', 'blue3'))
axis (1, seq (1.5, 9, 2.5),
      c ('Calves raised', 'Calves raised', 'Cows in', 'Calves raised'), tick= FALSE, line= 1.5)
axis (1, seq (1.5, 9, 2.5),
      c ('for red meat', 'for replacement', 'herds', 'for veal'), tick= FALSE, line= 2.5)
axis (1, 1:2, c ('beef', 'dairy'))
axis (1, 3.5:4.5, c ('beef', 'dairy'))
axis (1, 6:7, c ('beef', 'dairy'))
axis (1, 8.5:9.5, c ('beef', 'dairy'))
axis (2, labels= FALSE)

mtext ('Overall welfare risk', side= 2, las= 0, line= 2, font= 2)
mtext ('[scale of PC 1]', side= 2, las= 0, line= 1)
mtext ('lower\nwelfare\nrisk', side= 2, at= -5, line= 1)
mtext ('higher\nwelfare\nrisk', side= 2, at= 7, line= 1)
text (c (1, 2, 3.5, 4.5, 6, 7, 8.5, 9.5), rep (-6.2, 3),
      boxplot (Comp.1 ~ origin + categ, meat.df, plot= FALSE) [['n']], cex= 0.7)
text (0.5, -6.2, 'N=', cex= 0.7)

r <- 1
for (i in seq (1, 8.5, 2.5)) {
    lines (c (i, i+1), estim.C1.val [c (r, r+1), 'estim'], col= gray (0.2), lwd= 4)
    lines (c (i, i+1), estim.C1.val [c (r, r+1), 'lo.ci'], col= gray (0.2), lwd= 2)
    lines (c (i, i+1), estim.C1.val [c (r, r+1), 'up.ci'], col= gray (0.2), lwd= 2)
    r <- r + 2
}

dev.off ()


# Graphics other components
boxplot (Comp.2 ~ origin + categ, meat.df)
boxplot (Comp.3 ~ origin + categ, meat.df)
boxplot (Comp.4 ~ origin + categ, meat.df)

write.table (meat.df, file= 'meat.csv', quote= FALSE, sep= ';', dec= '.', row.names= FALSE)


## Mixed model extended version
#######################################################################################################

country.df <- read.table ('../Data/Country.csv', header= TRUE, sep= ';', dec= '.', as.is= FALSE)
row.names (country.df) <- country.df [, 'ResponseId']

country.df [country.df [, 'Q81'] == 'Austria' |
            country.df [, 'Q81'] == 'Czech Republic' |
            country.df [, 'Q81'] == 'Denmark' |
            country.df [, 'Q81'] == 'Finland' |
            country.df [, 'Q81'] == 'France' |
            country.df [, 'Q81'] == 'Germany' |
            country.df [, 'Q81'] == 'Ireland' |
            country.df [, 'Q81'] == 'Italy' |
            country.df [, 'Q81'] == 'Netherlands' |
            country.df [, 'Q81'] == 'Norway' |
            country.df [, 'Q81'] == 'Sweden' |
            country.df [, 'Q81'] == 'Switzerland' |
            country.df [, 'Q81'] == 'United Kingdom of Great Britain and Northern Ireland', 'region'] <- 'EU'
country.df [country.df [, 'Q81'] == 'Canada' |
            country.df [, 'Q81'] == 'United States of America', 'region'] <- 'NA'
country.df [country.df [, 'Q81'] == 'Australia' |
            country.df [, 'Q81'] == 'Brazil' |
            country.df [, 'Q81'] == 'Chile' |
            country.df [, 'Q81'] == 'Israel' |
            country.df [, 'Q81'] == 'Japan' |
            country.df [, 'Q81'] == 'Serbia' |
            country.df [, 'Q81'] == 'Turkey' |
            country.df [, 'Q81'] == 'Uruguay', 'region'] <- 'OT'
country.df [, 'region'] <- factor (country.df [, 'region'], levels= c ('EU', 'NA', 'OT'))
summary (country.df)

meat.df [, 'region'] <- country.df [meat.df [, 'id'], 'region']
summary (meat.df)
table (meat.df [, 'region'])

contrasts (meat.df [, 'origin']) <- contr.sum (2)
contrasts (meat.df [, 'categ']) <- contr.sum (4)
contrasts (meat.df [, 'region']) <- contr.sum (3)


meat.df <- cbind (meat.df, model.matrix (~ origin + categ + region, meat.df) [, -1])

meat.df [, 'oc1'] <- meat.df [, 'origin1'] * meat.df [, 'categ1']
meat.df [, 'oc2'] <- meat.df [, 'origin1'] * meat.df [, 'categ2']
meat.df [, 'oc3'] <- meat.df [, 'origin1'] * meat.df [, 'categ3']

meat.df [, 'or1'] <- meat.df [, 'origin1'] * meat.df [, 'region1']
meat.df [, 'or2'] <- meat.df [, 'origin1'] * meat.df [, 'region2']

meat.df [, 'cr1'] <- meat.df [, 'categ1'] * meat.df [, 'region1']
meat.df [, 'cr2'] <- meat.df [, 'categ1'] * meat.df [, 'region2']
meat.df [, 'cr3'] <- meat.df [, 'categ2'] * meat.df [, 'region1']
meat.df [, 'cr4'] <- meat.df [, 'categ2'] * meat.df [, 'region2']
meat.df [, 'cr5'] <- meat.df [, 'categ3'] * meat.df [, 'region1']
meat.df [, 'cr6'] <- meat.df [, 'categ3'] * meat.df [, 'region2']

meat.df [, 'ocr1'] <- meat.df [, 'origin1'] * meat.df [, 'categ1'] * meat.df [, 'region1']
meat.df [, 'ocr2'] <- meat.df [, 'origin1'] * meat.df [, 'categ1'] * meat.df [, 'region2']
meat.df [, 'ocr3'] <- meat.df [, 'origin1'] * meat.df [, 'categ2'] * meat.df [, 'region1']
meat.df [, 'ocr4'] <- meat.df [, 'origin1'] * meat.df [, 'categ2'] * meat.df [, 'region2']
meat.df [, 'ocr5'] <- meat.df [, 'origin1'] * meat.df [, 'categ3'] * meat.df [, 'region1']
meat.df [, 'ocr6'] <- meat.df [, 'origin1'] * meat.df [, 'categ3'] * meat.df [, 'region2']

names (meat.df)


# Maximum model
max.ext.lmer <- blmer (Comp.1 ~ ## Main effects
                       origin1 +
                       categ1 + categ2 + categ3 +
                       region1 + region2 +
                       ## two-way interactions
                       oc1 + oc2 + oc3 +
                       or1 + or2 +
                       cr1 + cr2 + cr3 + cr4 + cr5 + cr6 +
                       ## three-way interaction
#                       ocr1 + ocr2 + ocr3 + ocr4 + ocr5 + ocr6 +
                      (1 | id/categ),
                  meat.df, weights= confidence)
summary (max.ext.lmer)

# Residuals
max.simResid <- simulateResiduals (max.ext.lmer, 1000)
plot (max.simResid)
plotResiduals (max.simResid, meat.df [, 'origin'])
plotResiduals (max.simResid, meat.df [, 'categ'])
plotResiduals (max.simResid, meat.df [, 'region'])
qqnorm (resid (max.ext.lmer))


# p-values
int1.ext.lmer <- blmer (Comp.1 ~ ## Main effects
                       origin1 +
                       categ1 + categ2 + categ3 +
                       region1 + region2 +
                       ## two-way interactions
                       or1 + or2 +
                       cr1 + cr2 + cr3 + cr4 + cr5 + cr6 +
                       (1 | id/categ),
                   meat.df, weights= confidence)
#summary (int1.ext.lmer)
int1.p <- PBmodcomp (max.ext.lmer, int1.ext.lmer)
warnings ()
summary (int1.p)

int2.ext.lmer <- blmer (Comp.1 ~ ## Main effects
                       origin1 +
                       categ1 + categ2 + categ3 +
                       region1 + region2 +
                       ## two-way interactions
                       oc1 + oc2 + oc3 +
                       cr1 + cr2 + cr3 + cr4 + cr5 + cr6 +
                       (1 | id/categ),
                   meat.df, weights= confidence)
#summary (int2.ext.lmer)
int2.p <- PBmodcomp (max.ext.lmer, int2.ext.lmer)
warnings ()
summary (int2.p)

int3.ext.lmer <- blmer (Comp.1 ~ ## Main effects
                       origin1 +
                       categ1 + categ2 + categ3 +
                       region1 + region2 +
                       ## two-way interactions
                       oc1 + oc2 + oc3 +
                       or1 + or2 +
                       (1 | id/categ),
                   meat.df, weights= confidence)
#summary (int3.ext.lmer)
int3.p <- PBmodcomp (max.ext.lmer, int3.ext.lmer)
warnings ()
summary (int3.p)

BvD.ext.lmer <- blmer (Comp.1 ~ ## Main effects
                       categ1 + categ2 + categ3 +
                       region1 + region2 +
                       ## two-way interactions
                       oc1 + oc2 + oc3 +
                       or1 + or2 +
                       cr1 + cr2 + cr3 + cr4 + cr5 + cr6 +
                       (1 | id/categ),
                   meat.df, weights= confidence)
#summary (BvD.ext.lmer)
BvD.p <- PBmodcomp (max.ext.lmer, BvD.ext.lmer)
warnings ()
summary (BvD.p)

AC.ext.lmer <- blmer (Comp.1 ~ ## Main effects
                       origin1 +
                       region1 + region2 +
                       ## two-way interactions
                       oc1 + oc2 + oc3 +
                       or1 + or2 +
                       cr1 + cr2 + cr3 + cr4 + cr5 + cr6 +
                      (1 | id/categ),
                  meat.df, weights= confidence)
#summary (AC.ext.lmer)
AC.p <- PBmodcomp (max.ext.lmer, AC.ext.lmer)
warnings ()
summary (AC.p)

Reg.ext.lmer <- blmer (Comp.1 ~ ## Main effects
                       origin1 +
                       categ1 + categ2 + categ3 +
                       ## two-way interactions
                       oc1 + oc2 + oc3 +
                       or1 + or2 +
                       cr1 + cr2 + cr3 + cr4 + cr5 + cr6 +
                      (1 | id/categ),
                  meat.df, weights= confidence)
#summary (Reg.ext.lmer)
Reg.p <- PBmodcomp (max.ext.lmer, Reg.ext.lmer)
warnings ()
summary (Reg.p)


# Graphics Component 1
estim.ext.C1.df <- data.frame (origin1= rep (rep (c (1, -1), 4), rep (3, 8)),
                               categ1= rep (rep (c (1, 0, 0, -1), rep (2, 4)), rep (3, 8)),
                               categ2= rep (rep (c (0, 1, 0, -1), rep (2, 4)), rep (3, 8)),
                               categ3= rep (rep (c (0, 0, 1, -1), rep (2, 4)), rep (3, 8)),
                               region1= rep (c (1, 0, -1), 8),
                               region2= rep (c (0, 1, -1), 8))
estim.ext.C1.df [, 'oc1'] <- estim.ext.C1.df [, 'origin1'] * estim.ext.C1.df [, 'categ1']
estim.ext.C1.df [, 'oc2'] <- estim.ext.C1.df [, 'origin1'] * estim.ext.C1.df [, 'categ2']
estim.ext.C1.df [, 'oc3'] <- estim.ext.C1.df [, 'origin1'] * estim.ext.C1.df [, 'categ3']

estim.ext.C1.df [, 'or1'] <- estim.ext.C1.df [, 'origin1'] * estim.ext.C1.df [, 'region1']
estim.ext.C1.df [, 'or2'] <- estim.ext.C1.df [, 'origin1'] * estim.ext.C1.df [, 'region2']

estim.ext.C1.df [, 'cr1'] <- estim.ext.C1.df [, 'categ1'] * estim.ext.C1.df [, 'region1']
estim.ext.C1.df [, 'cr2'] <- estim.ext.C1.df [, 'categ1'] * estim.ext.C1.df [, 'region2']
estim.ext.C1.df [, 'cr3'] <- estim.ext.C1.df [, 'categ2'] * estim.ext.C1.df [, 'region1']
estim.ext.C1.df [, 'cr4'] <- estim.ext.C1.df [, 'categ2'] * estim.ext.C1.df [, 'region2']
estim.ext.C1.df [, 'cr5'] <- estim.ext.C1.df [, 'categ3'] * estim.ext.C1.df [, 'region1']
estim.ext.C1.df [, 'cr6'] <- estim.ext.C1.df [, 'categ3'] * estim.ext.C1.df [, 'region2']

estim.ext.C1.Mod <- function (x) predict (x, estim.ext.C1.df, re.form= NA)
estim.ext.C1.raw <- bootMer (max.ext.lmer, estim.ext.C1.Mod, nsim= 1000, .progress= 'win')
estim.ext.C1.val <- extract.ci (estim.ext.C1.raw)


postscript ('Mandel_C1_ext.eps', width= 16, height= 9, horizontal= FALSE, pointsize= 14, paper= 'special')
par (bty= 'n', las= 1, mar= c (5.5, 4, 0.1, 0.2))
boxplot (Comp.1 ~ region + origin + categ, meat.df,
         ylab= '', xlab= '', xaxt= 'n', yaxt= 'n',
         xlim= c (4/7, 25.5 + 3/7), at= c (1:6, seq (7.5, 12.5, 1), 14:19, seq (20.5, 25.5, 1)),
         range= 0, col= rep (c ('indianred', 'lightblue'), rep (3, 2)), whisklty= 1,
         whiskcol= rep (c ('red3', 'blue3'), rep (3, 2)),
         boxcol= rep (c ('red3', 'blue3'), rep (3, 2)),
         medcol= rep (c ('red3', 'blue3'), rep (3, 2)),
         staplecol= rep (c ('red3', 'blue3'), rep (3, 2)))

axis (1, seq (3.5, 23, 6.5),
      c ('Calves raised', 'Calves raised', 'Cows in', 'Calves raised'), tick= FALSE, line= 2.5)
axis (1, seq (3.5, 23, 6.5),
      c ('for red meat', 'for replacement', 'herds', 'for veal'), tick= FALSE, line= 3.5)

axis (1, c (2, 5), c ('beef', 'dairy'), tick= FALSE, line= 1.5)
axis (1, c (8.5, 11.5), c ('beef', 'dairy'), tick= FALSE, line= 1.5)
axis (1, c (15, 18), c ('beef', 'dairy'), tick= FALSE, line= 1.5)
axis (1, c (21.5, 24.5), c ('beef', 'dairy'), tick= FALSE, line= 1.5)

axis (1, 1:3, c ('EU', 'N-A', 'Other'), cex.axis= 0.8)
axis (1, 4:6, c ('EU', 'N-A', 'Other'), cex.axis= 0.8)
axis (1, seq (7.5, 9.5, 1), c ('EU', 'N-A', 'Other'), cex.axis= 0.8)
axis (1, seq (10.5, 12.5, 1), c ('EU', 'N-A', 'Other'), cex.axis= 0.8)
axis (1, 14:16, c ('EU', 'N-A', 'Other'), cex.axis= 0.8)
axis (1, 17:19, c ('EU', 'N-A', 'Other'), cex.axis= 0.8)
axis (1, seq (20.5, 22.5, 1), c ('EU', 'N-A', 'Other'), cex.axis= 0.8)
axis (1, seq (23.5, 25.5, 1), c ('EU', 'N-A', 'Other'), cex.axis= 0.8)

axis (2, labels= FALSE)

mtext ('Overall welfare risk', side= 2, las= 0, line= 2, font= 2)
mtext ('[scale of PC 1]', side= 2, las= 0, line= 1)
mtext ('lower\nwelfare\nrisk', side= 2, at= -5, line= 1)
mtext ('higher\nwelfare\nrisk', side= 2, at= 7, line= 1)
text (c (1:6, seq (7.5, 12.5, 1), 14:19, seq (20.5, 25.5, 1)),
      rep (-6.2, 24),
      boxplot (Comp.1 ~ region + origin + categ, meat.df, plot= FALSE) [['n']], cex= 0.7)
text (0.5, -6.2, 'N=', cex= 0.7)

r <- 1
for (i in seq (1, 20.5, 6.5)) {
    lines (seq (i, i+2, 1), estim.ext.C1.val [seq (r, r+2, 1), 'estim'], col= gray (0.2), lwd= 4)
    lines (seq (i, i+2, 1), estim.ext.C1.val [seq (r, r+2, 1), 'lo.ci'], col= gray (0.2), lwd= 2)
    lines (seq (i, i+2, 1), estim.ext.C1.val [seq (r, r+2, 1), 'up.ci'], col= gray (0.2), lwd= 2)
    r <- r + 6
}

r <- 4
for (i in seq (4, 23.5, 6.5)) {
    lines (seq (i, i+2, 1), estim.ext.C1.val [seq (r, r+2, 1), 'estim'], col= gray (0.2), lwd= 4)
    lines (seq (i, i+2, 1), estim.ext.C1.val [seq (r, r+2, 1), 'lo.ci'], col= gray (0.2), lwd= 2)
    lines (seq (i, i+2, 1), estim.ext.C1.val [seq (r, r+2, 1), 'up.ci'], col= gray (0.2), lwd= 2)
    r <- r + 6
}

dev.off ()
