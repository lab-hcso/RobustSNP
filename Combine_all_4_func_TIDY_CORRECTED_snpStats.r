require(car)
require(mvtnorm)
require(snpStats)
Rlinear <- function(data = NA, outcome = NA, SNP = NA,
    InputFile = NA, covPres = TRUE, COVAR = NA, CovarFile = NA,
    inputheader = TRUE, covarheader = TRUE, maxpts = 25000, abseps = 0.001) {
    if (is.na(SNP))
        SNP = 1:ncol(data)
    if (!is.na(InputFile)) {
        Rinput = read.table(InputFile, header = inputheader)
        outcome = Rinput[, 1]
        data = Rinput[, -1]
    }
    if (!is.na(CovarFile))
        COVAR = read.table(CovarFile, header = covarheader)
    n.subj = nrow(data)
    if (!covPres)
        x.adj <- rep(1, n.subj)
    else x.adj <- cbind(rep(1, n.subj), COVAR)
    x.adj = as.matrix(x.adj)
    reg.out <- glm.fit(y = outcome, x = x.adj, family = gaussian())
    mu <- reg.out$fitted.values
    resid = reg.out$residuals
    a <- sum(resid^2)/reg.out$df.residual
    v <- 1/a
    marker.score.max <- function(y, geno, maxpts = maxpts, abseps = abseps) {
        miss = which(is.na(geno))
        if (length(miss) > 0) {
            geno = geno[-miss]
            if (covPres)
                x.adj = x.adj[-miss, ]
            else x.adj = x.adj[-miss]
            y = y[-miss]
            mu = mu[-miss]
        }
        n.subj = length(geno)
        u.mtx <- (y - mu) * geno/a
        u.score <- sum(u.mtx)
        vscoreLin = sum(u.mtx * t(u.mtx))
        v.11 <- t(x.adj * v) %*% x.adj
        v.21 <- t(geno) %*% (x.adj * v)
        v.22 <- t(geno * v) %*% geno
        size = length(y)
        Sbeta = u.mtx
        Salpha.mtx = (y - mu) * x.adj/a
        eff.score.add = numeric(size)
        inv11 = solve(v.11)
        if (covPres) {
            for (i in 1:size) {
                eff.score.add[i] <- Sbeta[i] - v.21 %*% inv11 %*%
                  Salpha.mtx[i, ]
            }
        }
        if (!covPres) {
            for (i in 1:size) {
                eff.score.add[i] <- Sbeta[i] - v.21 %*% inv11 %*%
                  Salpha.mtx[i]
            }
        }
        sum.eff.score.add = sum(eff.score.add)
        VLin.add = sum(eff.score.add * t(eff.score.add))
        Z.add = sum.eff.score.add/sqrt(VLin.add)
        P.add = pnorm(abs(Z.add), lower = FALSE) * 2
        if ((sum(geno == 2) == 0)) {
            Z.rec = -9
            Z.dom = Z.add
            P.rec = -9
            P.dom = P.add
            theop = P.add
            integ.error = -9
            final.results = c(Z.add, Z.rec, Z.dom, P.add, P.rec,
                P.dom, theop, integ.error)
            return(final.results)
            break
        }
        if ((sum(geno == 0) == 0)) {
            Z.rec = Z.add
            Z.dom = -9
            P.rec = P.add
            P.dom = -9
            theop = P.add
            integ.error = -9
            final.results = c(Z.add, Z.rec, Z.dom, P.add, P.rec,
                P.dom, theop, integ.error)
            return(final.results)
            break
        }
        if ((sum(geno == 1) == 0)) {
            Z.rec = Z.add
            Z.dom = Z.add
            P.rec = P.add
            P.dom = P.add
            theop = P.add
            integ.error = -9
            final.results = c(Z.add, Z.rec, Z.dom, P.add, P.rec,
                P.dom, theop, integ.error)
            return(final.results)
            break
        }
        recsnp = recode(geno, "1=0;2=1")
        u.mtx.rec <- (y - mu) * recsnp/a
        u.score.rec <- sum(u.mtx.rec)
        vscoreLin = sum(u.mtx.rec * t(u.mtx.rec))
        v.11 <- t(x.adj * v) %*% x.adj
        v.21.rec <- t(recsnp) %*% (x.adj * v)
        v.22.rec <- t(recsnp * v) %*% recsnp
        size = length(y)
        Sbeta.rec = u.mtx.rec
        Salpha.mtx = (y - mu) * x.adj/a
        eff.score.rec = numeric(size)
        if (covPres) {
            for (i in 1:size) {
                eff.score.rec[i] <- Sbeta.rec[i] - v.21.rec %*%
                  inv11 %*% Salpha.mtx[i, ]
            }
        }
        if (!covPres) {
            for (i in 1:size) {
                eff.score.rec[i] <- Sbeta.rec[i] - v.21.rec %*%
                  inv11 %*% Salpha.mtx[i]
            }
        }
        sum.eff.score.rec = sum(eff.score.rec)
        VLin.rec = sum(eff.score.rec * t(eff.score.rec))
        Z.rec = sum.eff.score.rec/sqrt(VLin.rec)
        P.rec = pnorm(abs(Z.rec), lower = FALSE) * 2
        domsnp = recode(geno, "2=1")
        u.mtx.dom <- (y - mu) * domsnp/a
        u.score.dom <- sum(u.mtx.dom)
        vscoreLin = sum(u.mtx.dom * t(u.mtx.dom))
        v.11 <- t(x.adj * v) %*% x.adj
        v.21.dom <- t(domsnp) %*% (x.adj * v)
        v.22.dom <- t(domsnp * v) %*% domsnp
        size = length(y)
        Sbeta.dom = u.mtx.dom
        Salpha.mtx = (y - mu) * x.adj/a
        eff.score.dom = numeric(size)
        if (covPres) {
            for (i in 1:size) {
                eff.score.dom[i] <- Sbeta.dom[i] - v.21.dom %*%
                  inv11 %*% Salpha.mtx[i, ]
            }
        }
        if (!covPres) {
            for (i in 1:size) {
                eff.score.dom[i] <- Sbeta.dom[i] - v.21.dom %*%
                  inv11 %*% Salpha.mtx[i]
            }
        }
        sum.eff.score.dom = sum(eff.score.dom)
        VLin.dom = sum(eff.score.dom * t(eff.score.dom))
        Z.dom = sum.eff.score.dom/sqrt(VLin.dom)
        P.dom = pnorm(abs(Z.dom), lower = FALSE) * 2
        maxz = max(abs(Z.add), abs(Z.rec), abs(Z.dom))
        VLin.add.rec = sum(eff.score.add * t(eff.score.rec))
        VLin.add.dom = sum(eff.score.add * t(eff.score.dom))
        VLin.rec.dom = sum(eff.score.rec * t(eff.score.dom))
        Z.add.rec = VLin.add.rec/(sqrt(VLin.add * VLin.rec))
        Z.add.dom = VLin.add.dom/(sqrt(VLin.add * VLin.dom))
        Z.rec.dom = VLin.rec.dom/(sqrt(VLin.rec * VLin.dom))
        a = Z.add.rec
        b = Z.add.dom
        c = Z.rec.dom
        mean0 = c(0, 0, 0)
        cov.Z = matrix(c(1, a, b, a, 1, c, b, c, 1), nrow = 3,
            byrow = TRUE)
        theop.obj = pmvnorm(lower = rep(-maxz, 3), upper = rep(maxz,
            3), mean = mean0, sigma = cov.Z, maxpts = maxpts,
            abseps = abseps, releps = 0)
        theop = 1 - theop.obj[1]
        integ.error = attributes(theop.obj)$error
        final.results = c(Z.add, Z.rec, Z.dom, P.add, P.rec,
            P.dom, theop, integ.error)
        return(final.results)
    }
    nsnps = ncol(data)
    maxobj1 = matrix(nrow = nsnps, ncol = 8)
    for (i in 1:nsnps) {
        maxobj1[i, ] = marker.score.max(y = outcome, geno = data[,
            i], maxpts = maxpts, abseps = abseps)
    }
    colnames(maxobj1) <- c("Z.add", "Z.rec", "Z.dom", "P.add",
        "P.rec", "P.dom", "theoP", "integ.error")
    maxobj1 = as.data.frame(cbind(SNP, maxobj1))
    return(maxobj1)
}
require(car)
require(mvtnorm)
require(snpStats)
Rbin <- function(data = NA, outcome = NA, SNP = NA,
    InputFile = NA, covPres = TRUE, COVAR = NA, CovarFile = NA,
    inputheader = TRUE, covarheader = TRUE, maxpts = 25000, abseps = 0.001) {
    a = 1
    if (is.na(SNP))
        SNP = 1:ncol(data)
    if (!is.na(InputFile)) {
        Rinput = read.table(InputFile, header = inputheader)
        outcome = Rinput[, 1]
        data = Rinput[, -1]
    }
    if (!is.na(CovarFile))
        COVAR = read.table(CovarFile, header = covarheader)
    n.subj = nrow(data)
    if (!covPres)
        x.adj <- rep(1, n.subj)
    else x.adj <- cbind(rep(1, n.subj), COVAR)
    x.adj = as.matrix(x.adj)
    col1 <- rep(1, nrow(data))
    reg.out <- glm.fit(y = outcome, x = x.adj, family = binomial())
    mu <- reg.out$fitted.values
    v <- mu * (1 - mu)
    marker.score.max <- function(y, geno, trait.type = "binomial",
        maxpts = maxpts, abseps = abseps) {
        miss = which(is.na(geno))
        if (length(miss) > 0) {
            geno = geno[-miss]
            if (covPres)
                x.adj = x.adj[-miss, ]
            else x.adj = x.adj[-miss]
            y = y[-miss]
            mu = mu[-miss]
            v = v[-miss]
        }
        n.subj = length(geno)
        u.mtx <- (y - mu) * geno/a
        u.score <- sum(u.mtx)
        vscoreLin = sum(u.mtx * t(u.mtx))
        v.11 <- t(x.adj * v) %*% x.adj
        v.21 <- t(geno) %*% (x.adj * v)
        v.22 <- t(geno * v) %*% geno
        size = length(y)
        Sbeta = u.mtx
        Salpha.mtx = (y - mu) * x.adj/a
        eff.score.add = numeric(size)
        inv11 = solve(v.11)
        if (covPres) {
            for (i in 1:size) {
                eff.score.add[i] <- Sbeta[i] - v.21 %*% inv11 %*%
                  Salpha.mtx[i, ]
            }
        }
        if (!covPres) {
            for (i in 1:size) {
                eff.score.add[i] <- Sbeta[i] - v.21 %*% inv11 %*%
                  Salpha.mtx[i]
            }
        }
        sum.eff.score.add = sum(eff.score.add)
        VLin.add = sum(eff.score.add * t(eff.score.add))
        Z.add = sum.eff.score.add/sqrt(VLin.add)
        P.add = pnorm(abs(Z.add), lower = FALSE) * 2
        if ((sum(geno == 2) == 0)) {
            Z.rec = -9
            Z.dom = Z.add
            P.rec = -9
            P.dom = P.add
            theop = P.add
            integ.error = -9
            final.results = c(Z.add, Z.rec, Z.dom, P.add, P.rec,
                P.dom, theop, integ.error)
            return(final.results)
            break
        }
        if ((sum(geno == 0) == 0)) {
            Z.rec = Z.add
            Z.dom = -9
            P.rec = P.add
            P.dom = -9
            theop = P.add
            integ.error = -9
            final.results = c(Z.add, Z.rec, Z.dom, P.add, P.rec,
                P.dom, theop, integ.error)
            return(final.results)
            break
        }
        if ((sum(geno == 1) == 0)) {
            Z.rec = Z.add
            Z.dom = Z.add
            P.rec = P.add
            P.dom = P.add
            theop = P.add
            integ.error = -9
            final.results = c(Z.add, Z.rec, Z.dom, P.add, P.rec,
                P.dom, theop, integ.error)
            return(final.results)
            break
        }
        recsnp = recode(geno, "1=0;2=1")
        u.mtx.rec <- (y - mu) * recsnp/a
        u.score.rec <- sum(u.mtx.rec)
        vscoreLin = sum(u.mtx.rec * t(u.mtx.rec))
        v.11 <- t(x.adj * v) %*% x.adj
        v.21.rec <- t(recsnp) %*% (x.adj * v)
        v.22.rec <- t(recsnp * v) %*% recsnp
        size = length(y)
        Sbeta.rec = u.mtx.rec
        Salpha.mtx = (y - mu) * x.adj/a
        eff.score.rec = numeric(size)
        if (covPres) {
            for (i in 1:size) {
                eff.score.rec[i] <- Sbeta.rec[i] - v.21.rec %*%
                  inv11 %*% Salpha.mtx[i, ]
            }
        }
        if (!covPres) {
            for (i in 1:size) {
                eff.score.rec[i] <- Sbeta.rec[i] - v.21.rec %*%
                  inv11 %*% Salpha.mtx[i]
            }
        }
        sum.eff.score.rec = sum(eff.score.rec)
        VLin.rec = sum(eff.score.rec * t(eff.score.rec))
        Z.rec = sum.eff.score.rec/sqrt(VLin.rec)
        P.rec = pnorm(abs(Z.rec), lower = FALSE) * 2
        domsnp = recode(geno, "2=1")
        u.mtx.dom <- (y - mu) * domsnp/a
        u.score.dom <- sum(u.mtx.dom)
        vscoreLin = sum(u.mtx.dom * t(u.mtx.dom))
        v.11 <- t(x.adj * v) %*% x.adj
        v.21.dom <- t(domsnp) %*% (x.adj * v)
        v.22.dom <- t(domsnp * v) %*% domsnp
        size = length(y)
        Sbeta.dom = u.mtx.dom
        Salpha.mtx = (y - mu) * x.adj/a
        eff.score.dom = numeric(size)
        if (covPres) {
            for (i in 1:size) {
                eff.score.dom[i] <- Sbeta.dom[i] - v.21.dom %*%
                  inv11 %*% Salpha.mtx[i, ]
            }
        }
        if (!covPres) {
            for (i in 1:size) {
                eff.score.dom[i] <- Sbeta.dom[i] - v.21.dom %*%
                  inv11 %*% Salpha.mtx[i]
            }
        }
        sum.eff.score.dom = sum(eff.score.dom)
        VLin.dom = sum(eff.score.dom * t(eff.score.dom))
        Z.dom = sum.eff.score.dom/sqrt(VLin.dom)
        P.dom = pnorm(abs(Z.dom), lower = FALSE) * 2
        maxz = max(abs(Z.add), abs(Z.rec), abs(Z.dom))
        VLin.add.rec = sum(eff.score.add * t(eff.score.rec))
        VLin.add.dom = sum(eff.score.add * t(eff.score.dom))
        VLin.rec.dom = sum(eff.score.rec * t(eff.score.dom))
        Z.add.rec = VLin.add.rec/(sqrt(VLin.add * VLin.rec))
        Z.add.dom = VLin.add.dom/(sqrt(VLin.add * VLin.dom))
        Z.rec.dom = VLin.rec.dom/(sqrt(VLin.rec * VLin.dom))
        a = Z.add.rec
        b = Z.add.dom
        c = Z.rec.dom
        mean0 = c(0, 0, 0)
        cov.Z = matrix(c(1, a, b, a, 1, c, b, c, 1), nrow = 3,
            byrow = TRUE)
        theop.obj = pmvnorm(lower = rep(-maxz, 3), upper = rep(maxz,
            3), mean = mean0, sigma = cov.Z, maxpts = maxpts,
            abseps = abseps, releps = 0)
        theop = 1 - theop.obj[1]
        integ.error = attributes(theop.obj)$error
        final.results = c(Z.add, Z.rec, Z.dom, P.add, P.rec,
            P.dom, theop, integ.error)
        return(final.results)
    }
    nsnps = ncol(data)
    maxobj1 = matrix(nrow = nsnps, ncol = 8)
    for (i in 1:nsnps) {
        maxobj1[i, ] = marker.score.max(y = outcome, geno = data[,
            i], trait.type = "binomial", maxpts = maxpts, abseps = abseps)
    }
    colnames(maxobj1) <- c("Z.add", "Z.rec", "Z.dom", "P.add",
        "P.rec", "P.dom", "theoP", "integ.error")
    maxobj1 = as.data.frame(cbind(SNP, maxobj1))
    return(maxobj1)
}
require(car)
require(mvtnorm)
require(snpStats)
Rlinear.block <- function(GenoFile = NA, covPres = TRUE,
    COVAR = NA, CovarFile = NA, SNPrange = NA, covarheader = TRUE,
    blocksize = 5000, maxpts = 25000, abseps = 0.001) {
    plink.data = read.plink(GenoFile)
    if (!is.na(SNPrange[1]))
        plink.data = plink.data[, SNPrange]
    fam = read.table(paste(GenoFile, ".fam", sep = ""), header = FALSE)
    outcome = fam[, 6]
    bim = read.table(paste(GenoFile, ".bim", sep = ""), header = FALSE)
    if (!is.na(SNPrange[1]))
        SNP = as.character(bim$V2[SNPrange])
    else SNP = as.character(bim$V2)
    if (!is.na(CovarFile))
        COVAR = read.table(CovarFile, header = covarheader)
    n.subj = nrow(plink.data)
    if (!covPres)
        x.adj <- rep(1, n.subj)
    else x.adj <- cbind(rep(1, n.subj), COVAR)
    x.adj = as.matrix(x.adj)
    reg.out <- glm.fit(y = outcome, x = x.adj, family = gaussian())
    mu <- reg.out$fitted.values
    resid = reg.out$residuals
    a <- sum(resid^2)/reg.out$df.residual
    v <- 1/a
    marker.score.max <- function(y, geno, maxpts = maxpts, abseps = abseps) {
        miss = which(is.na(geno))
        if (length(miss) > 0) {
            geno = geno[-miss]
            if (covPres)
                x.adj = x.adj[-miss, ]
            else x.adj = x.adj[-miss]
            y = y[-miss]
            mu = mu[-miss]
        }
        n.subj = length(geno)
        u.mtx <- (y - mu) * geno/a
        u.score <- sum(u.mtx)
        vscoreLin = sum(u.mtx * t(u.mtx))
        v.11 <- t(x.adj * v) %*% x.adj
        v.21 <- t(geno) %*% (x.adj * v)
        v.22 <- t(geno * v) %*% geno
        size = length(y)
        Sbeta = u.mtx
        Salpha.mtx = (y - mu) * x.adj/a
        eff.score.add = numeric(size)
        inv11 = solve(v.11)
        if (covPres) {
            for (i in 1:size) {
                eff.score.add[i] <- Sbeta[i] - v.21 %*% inv11 %*%
                  Salpha.mtx[i, ]
            }
        }
        if (!covPres) {
            for (i in 1:size) {
                eff.score.add[i] <- Sbeta[i] - v.21 %*% inv11 %*%
                  Salpha.mtx[i]
            }
        }
        sum.eff.score.add = sum(eff.score.add)
        VLin.add = sum(eff.score.add * t(eff.score.add))
        Z.add = sum.eff.score.add/sqrt(VLin.add)
        P.add = pnorm(abs(Z.add), lower = FALSE) * 2
        if ((sum(geno == 2) == 0)) {
            Z.rec = -9
            Z.dom = Z.add
            P.rec = -9
            P.dom = P.add
            theop = P.add
            integ.error = -9
            final.results = c(Z.add, Z.rec, Z.dom, P.add, P.rec,
                P.dom, theop, integ.error)
            return(final.results)
            break
        }
        if ((sum(geno == 0) == 0)) {
            Z.rec = Z.add
            Z.dom = -9
            P.rec = P.add
            P.dom = -9
            theop = P.add
            integ.error = -9
            final.results = c(Z.add, Z.rec, Z.dom, P.add, P.rec,
                P.dom, theop, integ.error)
            return(final.results)
            break
        }
        if ((sum(geno == 1) == 0)) {
            Z.rec = Z.add
            Z.dom = Z.add
            P.rec = P.add
            P.dom = P.add
            theop = P.add
            integ.error = -9
            final.results = c(Z.add, Z.rec, Z.dom, P.add, P.rec,
                P.dom, theop, integ.error)
            return(final.results)
            break
        }
        recsnp = recode(geno, "1=0;2=1")
        u.mtx.rec <- (y - mu) * recsnp/a
        u.score.rec <- sum(u.mtx.rec)
        vscoreLin = sum(u.mtx.rec * t(u.mtx.rec))
        v.11 <- t(x.adj * v) %*% x.adj
        v.21.rec <- t(recsnp) %*% (x.adj * v)
        v.22.rec <- t(recsnp * v) %*% recsnp
        size = length(y)
        Sbeta.rec = u.mtx.rec
        Salpha.mtx = (y - mu) * x.adj/a
        eff.score.rec = numeric(size)
        if (covPres) {
            for (i in 1:size) {
                eff.score.rec[i] <- Sbeta.rec[i] - v.21.rec %*%
                  inv11 %*% Salpha.mtx[i, ]
            }
        }
        if (!covPres) {
            for (i in 1:size) {
                eff.score.rec[i] <- Sbeta.rec[i] - v.21.rec %*%
                  inv11 %*% Salpha.mtx[i]
            }
        }
        sum.eff.score.rec = sum(eff.score.rec)
        VLin.rec = sum(eff.score.rec * t(eff.score.rec))
        Z.rec = sum.eff.score.rec/sqrt(VLin.rec)
        P.rec = pnorm(abs(Z.rec), lower = FALSE) * 2
        domsnp = recode(geno, "2=1")
        u.mtx.dom <- (y - mu) * domsnp/a
        u.score.dom <- sum(u.mtx.dom)
        vscoreLin = sum(u.mtx.dom * t(u.mtx.dom))
        v.11 <- t(x.adj * v) %*% x.adj
        v.21.dom <- t(domsnp) %*% (x.adj * v)
        v.22.dom <- t(domsnp * v) %*% domsnp
        size = length(y)
        Sbeta.dom = u.mtx.dom
        Salpha.mtx = (y - mu) * x.adj/a
        eff.score.dom = numeric(size)
        if (covPres) {
            for (i in 1:size) {
                eff.score.dom[i] <- Sbeta.dom[i] - v.21.dom %*%
                  inv11 %*% Salpha.mtx[i, ]
            }
        }
        if (!covPres) {
            for (i in 1:size) {
                eff.score.dom[i] <- Sbeta.dom[i] - v.21.dom %*%
                  inv11 %*% Salpha.mtx[i]
            }
        }
        sum.eff.score.dom = sum(eff.score.dom)
        VLin.dom = sum(eff.score.dom * t(eff.score.dom))
        Z.dom = sum.eff.score.dom/sqrt(VLin.dom)
        P.dom = pnorm(abs(Z.dom), lower = FALSE) * 2
        maxz = max(abs(Z.add), abs(Z.rec), abs(Z.dom))
        VLin.add.rec = sum(eff.score.add * t(eff.score.rec))
        VLin.add.dom = sum(eff.score.add * t(eff.score.dom))
        VLin.rec.dom = sum(eff.score.rec * t(eff.score.dom))
        Z.add.rec = VLin.add.rec/(sqrt(VLin.add * VLin.rec))
        Z.add.dom = VLin.add.dom/(sqrt(VLin.add * VLin.dom))
        Z.rec.dom = VLin.rec.dom/(sqrt(VLin.rec * VLin.dom))
        a = Z.add.rec
        b = Z.add.dom
        c = Z.rec.dom
        mean0 = c(0, 0, 0)
        cov.Z = matrix(c(1, a, b, a, 1, c, b, c, 1), nrow = 3,
            byrow = TRUE)
        theop.obj = pmvnorm(lower = rep(-maxz, 3), upper = rep(maxz,
            3), mean = mean0, sigma = cov.Z, maxpts = maxpts,
            abseps = abseps, releps = 0)
        theop = 1 - theop.obj[1]
        integ.error = attributes(theop.obj)$error
        final.results = c(Z.add, Z.rec, Z.dom, P.add, P.rec,
            P.dom, theop, integ.error)
        return(final.results)
    }
    nsnps = ncol(plink.data)
    maxobj1 = matrix(nrow = nsnps, ncol = 8)
    no.block = floor(nsnps/blocksize) + 1
    endSNPpart1 = 1:(no.block - 1) * blocksize
    endSNP = c(1:(no.block - 1) * blocksize, nsnps)
    startSNP = c(1, endSNPpart1 + 1)
    for (B in 1:no.block) {
        data = as(plink.data[, startSNP[B]:endSNP[B]], "numeric")
        for (i in startSNP[B]:endSNP[B]) {
            maxobj1[i, ] = marker.score.max(y = outcome, geno = data[,
                i - startSNP[B] + 1], trait.type = "binomial",
                maxpts = maxpts, abseps = abseps)
        }
    }
    colnames(maxobj1) <- c("Z.add", "Z.rec", "Z.dom", "P.add",
        "P.rec", "P.dom", "theoP", "integ.error")
    maxobj1 = as.data.frame(cbind(SNP, maxobj1))
    return(maxobj1)
}
require(car)
require(mvtnorm)
require(snpStats)
Rbin.block <- function(GenoFile = NA, covPres = TRUE,
    COVAR = NA, CovarFile = NA, SNPrange = NA, covarheader = TRUE,
    blocksize = 5000, maxpts = 25000, abseps = 0.001) {
    a = 1
    plink.data = read.plink(GenoFile)
    if (!is.na(SNPrange[1]))
        plink.data = plink.data[, SNPrange]
    fam = read.table(paste(GenoFile, ".fam", sep = ""), header = FALSE)
    outcome = fam[, 6] - 1
    bim = read.table(paste(GenoFile, ".bim", sep = ""), header = FALSE)
    if (!is.na(SNPrange[1]))
        SNP = as.character(bim$V2[SNPrange])
    else SNP = as.character(bim$V2)
    if (!is.na(CovarFile))
        COVAR = read.table(CovarFile, header = covarheader)
    n.subj = nrow(plink.data)
    if (!covPres)
        x.adj <- rep(1, n.subj)
    else x.adj <- cbind(rep(1, n.subj), COVAR)
    x.adj = as.matrix(x.adj)
    reg.out <- glm.fit(y = outcome, x = x.adj, family = binomial())
    mu <- reg.out$fitted.values
    v <- mu * (1 - mu)
    marker.score.max <- function(y, geno, trait.type = "binomial",
        maxpts = maxpts, abseps = abseps) {
        miss = which(is.na(geno))
        if (length(miss) > 0) {
            geno = geno[-miss]
            if (!is.na(CovarFile))
                x.adj = x.adj[-miss, ]
            else x.adj = x.adj[-miss]
            y = y[-miss]
            mu = mu[-miss]
            v = v[-miss]
        }
        n.subj = length(geno)
        u.mtx <- (y - mu) * geno/a
        u.score <- sum(u.mtx)
        vscoreLin = sum(u.mtx * t(u.mtx))
        v.11 <- t(x.adj * v) %*% x.adj
        v.21 <- t(geno) %*% (x.adj * v)
        v.22 <- t(geno * v) %*% geno
        size = length(y)
        Sbeta = u.mtx
        Salpha.mtx = (y - mu) * x.adj/a
        eff.score.add = numeric(size)
        inv11 = solve(v.11)
        if (covPres) {
            for (i in 1:size) {
                eff.score.add[i] <- Sbeta[i] - v.21 %*% inv11 %*%
                  Salpha.mtx[i, ]
            }
        }
        if (!covPres) {
            for (i in 1:size) {
                eff.score.add[i] <- Sbeta[i] - v.21 %*% inv11 %*%
                  Salpha.mtx[i]
            }
        }
        sum.eff.score.add = sum(eff.score.add)
        VLin.add = sum(eff.score.add * t(eff.score.add))
        Z.add = sum.eff.score.add/sqrt(VLin.add)
        P.add = pnorm(abs(Z.add), lower = FALSE) * 2
        if ((sum(geno == 2) == 0)) {
            Z.rec = -9
            Z.dom = Z.add
            P.rec = -9
            P.dom = P.add
            theop = P.add
            integ.error = -9
            final.results = c(Z.add, Z.rec, Z.dom, P.add, P.rec,
                P.dom, theop, integ.error)
            return(final.results)
            break
        }
        if ((sum(geno == 0) == 0)) {
            Z.rec = Z.add
            Z.dom = -9
            P.rec = P.add
            P.dom = -9
            theop = P.add
            integ.error = -9
            final.results = c(Z.add, Z.rec, Z.dom, P.add, P.rec,
                P.dom, theop, integ.error)
            return(final.results)
            break
        }
        if ((sum(geno == 1) == 0)) {
            Z.rec = Z.add
            Z.dom = Z.add
            P.rec = P.add
            P.dom = P.add
            theop = P.add
            integ.error = -9
            final.results = c(Z.add, Z.rec, Z.dom, P.add, P.rec,
                P.dom, theop, integ.error)
            return(final.results)
            break
        }
        recsnp = recode(geno, "1=0;2=1")
        u.mtx.rec <- (y - mu) * recsnp/a
        u.score.rec <- sum(u.mtx.rec)
        vscoreLin = sum(u.mtx.rec * t(u.mtx.rec))
        v.11 <- t(x.adj * v) %*% x.adj
        v.21.rec <- t(recsnp) %*% (x.adj * v)
        v.22.rec <- t(recsnp * v) %*% recsnp
        size = length(y)
        Sbeta.rec = u.mtx.rec
        Salpha.mtx = (y - mu) * x.adj/a
        eff.score.rec = numeric(size)
        if (covPres) {
            for (i in 1:size) {
                eff.score.rec[i] <- Sbeta.rec[i] - v.21.rec %*%
                  inv11 %*% Salpha.mtx[i, ]
            }
        }
        if (!covPres) {
            for (i in 1:size) {
                eff.score.rec[i] <- Sbeta.rec[i] - v.21.rec %*%
                  inv11 %*% Salpha.mtx[i]
            }
        }
        sum.eff.score.rec = sum(eff.score.rec)
        VLin.rec = sum(eff.score.rec * t(eff.score.rec))
        Z.rec = sum.eff.score.rec/sqrt(VLin.rec)
        P.rec = pnorm(abs(Z.rec), lower = FALSE) * 2
        domsnp = recode(geno, "2=1")
        u.mtx.dom <- (y - mu) * domsnp/a
        u.score.dom <- sum(u.mtx.dom)
        vscoreLin = sum(u.mtx.dom * t(u.mtx.dom))
        v.11 <- t(x.adj * v) %*% x.adj
        v.21.dom <- t(domsnp) %*% (x.adj * v)
        v.22.dom <- t(domsnp * v) %*% domsnp
        size = length(y)
        Sbeta.dom = u.mtx.dom
        Salpha.mtx = (y - mu) * x.adj/a
        eff.score.dom = numeric(size)
        if (covPres) {
            for (i in 1:size) {
                eff.score.dom[i] <- Sbeta.dom[i] - v.21.dom %*%
                  inv11 %*% Salpha.mtx[i, ]
            }
        }
        if (!covPres) {
            for (i in 1:size) {
                eff.score.dom[i] <- Sbeta.dom[i] - v.21.dom %*%
                  inv11 %*% Salpha.mtx[i]
            }
        }
        sum.eff.score.dom = sum(eff.score.dom)
        VLin.dom = sum(eff.score.dom * t(eff.score.dom))
        Z.dom = sum.eff.score.dom/sqrt(VLin.dom)
        P.dom = pnorm(abs(Z.dom), lower = FALSE) * 2
        maxz = max(abs(Z.add), abs(Z.rec), abs(Z.dom))
        VLin.add.rec = sum(eff.score.add * t(eff.score.rec))
        VLin.add.dom = sum(eff.score.add * t(eff.score.dom))
        VLin.rec.dom = sum(eff.score.rec * t(eff.score.dom))
        Z.add.rec = VLin.add.rec/(sqrt(VLin.add * VLin.rec))
        Z.add.dom = VLin.add.dom/(sqrt(VLin.add * VLin.dom))
        Z.rec.dom = VLin.rec.dom/(sqrt(VLin.rec * VLin.dom))
        a = Z.add.rec
        b = Z.add.dom
        c = Z.rec.dom
        mean0 = c(0, 0, 0)
        cov.Z = matrix(c(1, a, b, a, 1, c, b, c, 1), nrow = 3,
            byrow = TRUE)
        theop.obj = pmvnorm(lower = rep(-maxz, 3), upper = rep(maxz,
            3), mean = mean0, sigma = cov.Z, maxpts = maxpts,
            abseps = abseps, releps = 0)
        theop = 1 - theop.obj[1]
        integ.error = attributes(theop.obj)$error
        final.results = c(Z.add, Z.rec, Z.dom, P.add, P.rec,
            P.dom, theop, integ.error)
        return(final.results)
    }
    nsnps = ncol(plink.data)
    maxobj1 = matrix(nrow = nsnps, ncol = 8)
    no.block = floor(nsnps/blocksize) + 1
    endSNPpart1 = 1:(no.block - 1) * blocksize
    endSNP = c(1:(no.block - 1) * blocksize, nsnps)
    startSNP = c(1, endSNPpart1 + 1)
    for (B in 1:no.block) {
        data = as(plink.data[, startSNP[B]:endSNP[B]], "numeric")
        for (i in startSNP[B]:endSNP[B]) {
            maxobj1[i, ] = marker.score.max(y = outcome, geno = data[,
                i - startSNP[B] + 1], trait.type = "binomial",
                maxpts = maxpts, abseps = abseps)
        }
    }
    colnames(maxobj1) <- c("Z.add", "Z.rec", "Z.dom", "P.add",
        "P.rec", "P.dom", "theoP", "integ.error")
    maxobj1 = as.data.frame(cbind(SNP, maxobj1))
    return(maxobj1)
}
