##
rholist <- c(0.6, 0.7, 0.75, 0.54, 0.45, 0.5, 0.01, 0.6, 0.8, 0.2, 
             0.4, 0.5, 0.6, 0.4, 0.8, 0.2, 0.1)

distlist <- c("norm", "norm", "norm", "pois", "pois", "pois", 
              "norm", "negbinom", "binom", "pois", "pois", 
              "binom", "negbinom", "negbinom", "pois", "norm", "norm")

nameslist <- c("attendanceTotal", "assessmentMath", "assessmentRead", 
               "attendance30day", "courseGradesCore", "courseGradesAll", 
               "GPA", "retention", "remedialCourse", "tardies", 
               "majorDiscipline", "minorDiscipline", "officeReferral", 
               "expulsionDays", "suspensionDays", "schoolMoves", "districtMoves")

struc <- list(dist = distlist, 
              rho = rholist,
              names = nameslist)

simplestruc <- list(dist = c("norm", "pois", "binom"), 
                    rho = c(0.2, 0.5, -0.2), 
                    names = c("test1", "test2", "test3"))

simplestruc <- list(dist = distlist[1:17], 
                    rho = rholist[1:17], 
                    names = nameslist[1:17])


dat <- genNumeric(1000, pattern = simplestruc, na.rm = TRUE, scale = .5)
summary(dat)
head(dat)

nvars <- ncol(dat)
modterms <- list(coefs = runif(nvars, min=-1, max=1), vars = sample(names(dat), nvars)) 

y <- genBinomialDV(dat, form=modterms, intercept=10)
dat <- cbind(y, dat)
table(dat$y)

