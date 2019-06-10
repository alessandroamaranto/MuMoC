load_flow_s <- function(shape_path, pCode, start, end) {
  
  library(dataRetrieval)
  library(dplyr)
  library(rgdal)
  
  crswgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  hp <- readOGR(shape_path)
  hp <- spTransform(hp, CRS = crswgs84)
  states <- c('NE', 'WY', 'TX', 'KS', 
              'SD', 'OK', 'CO', 'NM')
  
  Qs <- list()
  Ns <- length(states)
  for (i in 1:Ns) {
    Qs[[i]] <-
      readNWISdata(stateCd = states[i],
                   parameterCd = pCode,
                   service = 'site',
                   seriesCatalogOutput = TRUE,
                   siteStatus = 'active')
  }
  
  Qs <- do.call(rbind, Qs)
  Qs <- filter_q(Qs, start, end, pCode)
  
  coordinates(Qs) <- ~dec_long_va + dec_lat_va
  proj4string(Qs) <- proj4string(hp)
  
  idx <- over(Qs, hp)  
  Qs <- Qs[which(!is.na(idx)), c(2, 20, 22)]  
  
  return(Qs)
}

filter_q <- function(x, start, end, pCode) {
  
  filter(x, parm_cd %in% pCode) %>%
    filter(count_nu > 10000) %>%
    filter(as.Date(begin_date) <= as.Date(start)) %>%
    filter(as.Date(end_date) >= as.Date(end)) %>%
    mutate(period = as.Date(end_date) - as.Date(begin_date)) %>%
    filter(period > 15 * 365)
  
}

load_flow <- function(S, p, start, end) {
  
  S <- S@data
  N <- nrow(S)
  X <- list()
  for (i in 1:N) {
    X[[i]] <- readNWISdv(siteNumbers = S[i ,1],
                         parameterCd = p,
                         startDate = start,
                         endDate = end)
  }
  
  X <- Filter(Negate(is.null), X)
  X <- lapply(X, monthly_agg, start, end)
  
  r <- lapply(X, station_stats)
  r <- do.call(rbind, r)
  r <- which(r[ ,2] == 0)
  X <- X[r]
  
  return(X)
} 

monthly_agg <- function(x, start, end) {
  
  start <- as.Date(start)
  end <- as.Date(end)
  
  x <- x[, -c(1, 5)]
  Y <- format(as.Date(x[, 2]),
              "%Y")
  M <- format(as.Date(x[, 2]),
              "%m")
  
  x <- cbind(x, Y, M)
  colnames(x)[3] <- 'flow'
  
  y <- aggregate(flow ~ M + Y,
                 x,
                 function(z)
                   mean(z, na.rm = T))
  
  Date <- as.Date(paste(y[, 2],
                        y[, 1],
                        "01",
                        sep = "-"))
  
  y <- cbind(Date, y[, 3])
  Tp <- data.frame(Date = seq.Date(start, end, "months"))
  
  y <- merge(y, Tp, by = "Date", all = T)
  y <- data.frame(ID = rep(x[1, 1], nrow(y)),
                  Date = y[, 1],
                  flow = y[, 2])
  return(y)
}

station_stats <- function(x) {
  
  N <- nrow(x)
  L <- sum(is.na(x[ ,3])) 
  R <- L/N
  s <- data.frame(L = L, R = R)
  
  return(s)
}

exclude_meta <- function(m, x) {
  
  N <- length(x)
  n <- character(N)
  
  for (i in 1:N) {
    n[i] <- paste(x[[i]][1, 1])
  }
  
  idx <- which(m@data[ ,1] %in% n)
  m <- m[idx, ]
  return(m)
}

min_distance <- function(x, y) {
  
  library(geosphere)
  N <- nrow(x)
  d <- matrix(NA, N, 2)
  
  for (i in 1:N) {
    
    dis <- distCosine(x[i, ], y)
    d[i, 1] <- which.min(dis)
    d[i, 2] <- min(dis)
    
  }
  return(d)
}

river_coupling <- function(x, y, e){
  
  if (e == F) {
    Nh <- 4
    z <- couple_dates(x, y, Nh)
    No <- nrow(z)
    r <- list()
    
    # Iterate along time steps
    for (j in 1 : (Nh + 1)) {
      r[[j]] <- autoregressive_river_terms(j, z, Nh, No)
    }
    
    r <- do.call(cbind, r)
    r <- r[ ,c(1, 6:2)]
    colnames(r)[2:6] <- paste0("Qt", 4:0)
    
    idx <- which(r[ ,1] %in% x[ ,1])  
    r <- r[idx, ]  
  }else{
    r <- NULL
  }
  return(r)
}

couple_dates <- function(x, y, Nh) {
  
  N <- nrow(x)
  start <- x[1, 1]
  end <- x[N, 1]
  
  y[, 2] <- num2date(y[ ,2])
  
  Rstart <- which(y[ ,2] == start)
  Rend <- which(y[ ,2] == end)
  
  z <- y[(Rstart - Nh):Rend, ]
  
  return(z)
}

num2date <- function(x) {
  
  x <-  as.Date('1970-01-01') + x
  return(x)
  
}

autoregressive_river_terms <- function(j, z, Nh, No) {
  
  if (j == 1) {
    r <- z[-c(1:(Nh)) ,c(2, 3)]
  }else if (j == Nh + 1) {
    r <- data.frame(Qt4 = z[-c((No - Nh + 1) : No), 3])
  } else {
    r <- data.frame( Qt =  z[-c( 1:(Nh - j + 1),(No - j + 2) : No), 3])
  }
  return(r)
}

add_river <- function(x, s , e) {
  
  e <- as.character(e)
  switch(e,
         'FALSE' = {x <- normalize(x, s);
         }, 'TRUE' = {warning('river not found');},
         stop("specify case"))
  return(x)
  
}

normalize <- function(x, idx) {
  
  p1 <- 0.1; p2 <- 0.9
  x1 <- x[idx, ]
  x2 <- x[-idx, ]
  
  lb <- apply(x1, 2, min)
  ub <- apply(x1, 2, max)
  Nv <- ncol(x1)
  
  for (i in 1:Nv) {
    x1[ ,i] <- (((p2 - p1) *  (x1[ ,i] - lb[i]))/ (ub[i] - lb[i])) + p1
    x2[ ,i] <- (((p2 - p1) *  (x2[ ,i] - lb[i]))/ (ub[i] - lb[i])) + p1
  }
  
  return(list(x1, x2))
}

attach_date <- function(X, x, idx, case) {
  
  switch(case,
         'calibration'={x <- cbind(X[idx, 1], x);
         }, 'validation' = {x <- cbind(X[-idx, 1], x);},
         stop("specify case"))
  return(x)
  
}

select_columns <- function(ts, x, c) {
  
  Nvar <- ncol(x)
  y <- x[ , c(c - ts + 1, Nvar)]
  return(y)
  
}

input_variable_selection <- function(cand, xc, xv, e, Fl, rc, rv, j, selection) {
  
  if (selection == 'run') {
    
    
    library(foreach)
    library(doParallel)
    warning('computationally expensive, consider using norun')
    warnings()
    
    registerDoParallel(3)
    CI <- length(cand)
    Nh <- 4
    eps <- matrix(0, CI, Nh)
    
    for (u in 1:CI) {
      
      Nvar <- ncol(xc)
      id <- as.numeric(cand[[u]])
      xci <- xc[ ,c(id, (Nvar-3):Nvar)]
      xvi <- xv[ ,c(id, (Nvar-3):Nvar)]
      
      e <- as.character(e)
      
      switch(e,
             'FALSE' = {xci <- cbind(rc[ ,ncol(rc)-j], xci);
             }, 'TRUE' = {warning('river not found');},
             stop("specify case"))
      switch(e,
             'FALSE' = {xvi <- cbind(rv[ ,ncol(rv)-j], xvi);
             }, 'TRUE' = {warning('river not found');},
             stop("specify case"))
      
      
      cv_error <- function(x_c, f, m, ns) {
        
        library(RSNNS)
        library(hydroGOF)
        
        Nvar <- ncol(x_c)
        x_c <- as.data.frame(x_c)
        Y <- x_c[ ,Nvar]
        X <- x_c[,-Nvar]
        
        cf <- f[[m]]
        
        YY <- Y[cf]
        XX <- X[cf, ]
        
        #New Calibration
        Y <- Y[-cf]
        X <- X[-cf, ]
        
        CandidatesStruct<- 7
        
        itmax <- 500
        set.seed(42)
        ANN <- mlp(X, Y, maxit = 500, size = 7,
                   initFunc = "Randomize_Weights",
                   hiddenActFunc="Act_Logistic",
                   learnFunc = "Rprop",
                   outputActFunc = "Act_Logistic",
                   inputsTest = XX,
                   targetsTest = YY)
        
        
        out <- predict(ANN, XX)
        err <- NSE(out[, 1], YY)
        return(err)
      }
      
      ptime = system.time ({
        eps_u <- foreach(m = 1:10) %dopar%  cv_error(xci, Fl, m, -100)
      })
      eps[u, j] <- mean(do.call(rbind, eps_u))
    }
    
    nv <- ncol(xci)
    x_ar <- xci[ ,(nv-3):nv]
    
    ptime = system.time ({
      eps_ar <- foreach(m = 1:10) %dopar%  cv_error(x_ar, Fl, m, -100)
    })
    
    eps_ar <- mean(unlist(eps_ar))
    
    iv_idx <- which.max(eps[ ,j])
    dns <- eps[iv_idx ,j] - eps_ar
  }else {
    load('Data_in/IVS_out.RData')
  }
  return(list(iv_idx, dns))
}

ann_err <- function(cand, xc, xv, e, Fl, rc, rv, j, Vs) {
  
  id <- as.numeric(cand[[Vs]])
  Nvar <- ncol(xc)
  
  xc <- xc[ ,c(id, (Nvar-3):Nvar)]
  xv <- xv[ ,c(id, (Nvar-3):Nvar)]
  
  e <- as.character(e)
  
  switch(e,
         'FALSE' = {xc <- cbind(rc[ ,ncol(rc)-j], xc);
         }, 'TRUE' = {warning('river not found');},
         stop("specify case"))
  switch(e,
         'FALSE' = {xv <- cbind(rv[ ,ncol(rv)-j], xv);
         }, 'TRUE' = {warning('river not found');},
         stop("specify case"))
  
  ann_return <- function(x_c, f, m, ns) {
    
    library(RSNNS)
    library(hydroGOF)
    
    Nvar <- ncol(x_c)
    x_c <- as.data.frame(x_c)
    Y <- x_c[ ,Nvar]
    X <- x_c[,-Nvar]
    
    cf <- f[[m]]
    
    YY <- Y[cf]
    XX <- X[cf, ]
    
    #New Calibration
    Y <- Y[-cf]
    X <- X[-cf, ]
    
    Nnodes <- seq(3, 13, 2)
    it <- seq(300, 700, 200)
    
    n1 <- length(Nnodes)
    n2 <- length(it)
    dummy_nse <- -Inf
    
    for (p in 1:n1) {
      for (w in 1:n2) {
        set.seed(42)
        ANN <- mlp(X, Y, maxit = Nnodes[p], 
                   size = it[w],
                   initFunc = "Randomize_Weights",
                   hiddenActFunc="Act_Logistic",
                   learnFunc = "Rprop",
                   outputActFunc = "Act_Logistic",
                   inputsTest = XX,
                   targetsTest = YY)
        
        out <- predict(ANN, XX)
        err <- NSE(out[, 1], YY)
        
        if (err > dummy_nse) {
          model <- ANN
          dummy_nse <- err
        }
      }
    }
    return(model)
    
  }
  
  ptime = system.time ({
    ann <- foreach(m = 1:10) %dopar%  ann_return(train, folds, m, -100)
  })
  
  E <- ens_average(ann, xv)
  return(E)
}

ens_average <- function(model, xv) {
  
  Nens <- length(ann)
  Nt <- nrow(xv)
  Nvar <- ncol(xv)
  Y <- matrix(0, Nt, Nens)
  
  for (t in 1:Nens) {
    Y[ ,t] <- predict(model[[t]], xv[ ,-Nvar])
  }
  
  Y <- apply(Y, 1, mean)
  error <- NSE(Y, xv[ ,Nvar])
  return(error)
}

table_3 <- function(x, s, shape_path) {
  
  library(rgdal)
  crswgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  
  hp.states <- c("South Dakota",
                 "Oklahoma",
                 "New Mexico",
                 "Nebraska",
                 "Texas",
                 "Wyoming",
                 "Colorado",
                 "Kansas")
  
  us <- readOGR(shape_path)
  us <- us[which(us@data$STATE_NAME %in% hp.states), ]
  us <- spTransform(us, crswgs84)
  coordinates(s) <- ~dec_long_va + dec_lat_va
  s@proj4string <- crswgs84
  
  Ns <- length(hp.states); Nt <- ncol(x)
  e_s <- matrix(0, Ns, Nt+1)
  
  for (i in 1:Ns) {
    
    idx <- which(us@data$STATE_NAME %in% hp.states[i])
    state <- us[idx, ]
    idx_s <- which(!is.na (over(s, state,  fn = NULL))[ ,1])
    
    x_s <- x[idx_s, ] 
    
    if (length(x_s) > 4 ) {
      e_s[i, -1] <- round(apply(x_s, 2, mean), 2)
      
    }else{
      e_s[i, 2:5] <- round(x_s, 2)
      
    }
    e_s[i, 1] <- length(idx_s)
  }
  return(e_s)
}

exhaustive_pareto_search <- function(b, X, s, e_v) {
  
  Nh <- ncol(e_v)
  Pf <- list()
  dm <- c(100, 50) # max distance a neighbour is found (save ct)
  
  
  for (k in 1:Nh) {
    print(k)
    # Select bad models
    bh <- b[[k]]
    sb <- s[bh, 5:4]; Nb <- nrow(sb)
    
    # Allocate memory for pareto front for specific lead times
    Pw <- list()
    
    # Iterate along them
    for (i in 1:Nb) {
      print(i)
      # Extract coordinates of bad well
      p <- sb[i, ]
      d <- distCosine(p, s[,5:4], r = 6378137)/1000
      d[bh] <- NA
      
      if (p[2] > 39) {
        idx <- which(d < dm[1])
      }else{
        idx <- which(d < dm[2])
      }
      
      if (length(idx) == 0 ) {
        d.small = min(d, na.rm = T)
        idx <- which(d < (d.small) + d.small*0.15)
      }
      
      NN <- length(idx)
      
      J1 <- e_v[idx, k]
      J2 <- vector(length = NN)
      J3 <- vector(length = NN)
      
      x <- X[[bh[i]]][ ,c(1, 43)]
      
      # Iterate along neighbouring wells
      for (j in 1:NN) {
        
        # Take data from neighbouring wells
        xnj <- X[[idx[j]]][, c(1, 43)]
        xnj <- merge(x, xnj,
                     by = "Date", all = T)
        
        # Extract only common dates
        idx2 <- which(xnj$Date %in% x$Date)
        xnj <- xnj[idx2, ]
        
        # Compute Obj2
        J2[j] <- abs(ccf(xn[ ,2],
                         xnj[ ,3],
                         na.action = na.pass,
                         lag.max = 0,
                         plot = F)$acf)
        
        J3[j] <- length(which(is.na(xnj$Wlt_0.y)) )
      }
      
      # Removes row with NA
      J2[which(is.na(J2))] = 0
      J <- data.frame(NSE = J1, Corr = J2, Na = J3, Wid = idx)
      rm_idx <- which(J[ ,3] == 0)
      J <- J[-rm_idx, ]
      
      # Add constraints
      Jc <- J
      Cr <- which(Jc[ ,2] < 0.5)
      Jc[Cr, ] <- 0
      Ceps <- which(Jc[ ,1] < 0.7)
      Jc[Ceps, ] <- 0
      
      # Compute pareto front
      P <- psel(Jc, high(NSE) * high(Corr))
      if (nrow(P) < 5) {
        P <- psel(Jc, high(NSE) * high(Corr),
                  top = min(5, length(which(Jc[ ,4] != 0)))
                  )
      }
      Pw[[i]] <- P
    }
    # Save all pareto fronts in a big list
    Pf[[k]] <- Pw
  }
  return(Pf)
}

NeighTrainTest <- function(x, y, z, flag) {
  
  
  # Split train and test set for neighbouring wells
  common.dates <- which(x[ ,1] %in%  y)
  train <- x[common.dates, ]
  test.temp <- x[-common.dates, ]
  
  # Perform the same operation with original data
  common.bad <- which(z[ ,1] %in% y)
  train.bad <- z[common.bad, ]
  test.bad <- z[-common.bad, ]
  
  # Compute test set for neighbouring well-
  test <-  test.temp[which(test.temp[ ,1] %in% test.bad[ ,1]), ]
  
  if (flag == 1) {
    return(train)
  }else {
    return(test)
  }
}

MCReturn <- function(x, y, flag) {
  
  if (flag == 1) {
    
    # Import libraries
    library(RWeka)
    
    # Rename columns and run the model
    colnames(x)[ncol(x)] <- "out"
    SelModel <- M5P(out~., data = x)
    
    # knn
  }else if(flag == 2){
    
    # Select max number of neighbours
    n.neigh <- 30
    
    # Prepare the data: input
    x.tr <- x[ ,-ncol(x)]
    x.te <- y[ ,-ncol(x)]
    
    # output
    y.tr <- x[ ,ncol(x)]
    
    # Fake NSE
    ns <- -1000
    
    for (g in 1:n.neigh) {
      k.nn <- knn.reg(train = x.tr, test = x.te, y = y.tr, k = g, algorithm = "kd_tree")
      NS <- NSE(k.nn$pred, y[ ,ncol(y)])
      if (NS > ns) {
        ns <- NS
        SelModel <- k.nn
      }
    }
    #Lwr
  }else if(flag == 3){
    
    cands.lwr <- seq(0.1, 0.9, 0.1)
    deg <- c(0, 1, 2)
    
    colnames(x)[ncol(x)] <- "out"
    # Check number of predictors
    
    # Fake NSE
    ns <- -1000
    
    if (ncol(x) <= 5) {
      for (g in 1:length(cands.lwr)) {
        for (p in 1:length(deg)) {
          lwr <- loess(out~., data = x,  model = T,
                       span = cands.lwr[g], degree = deg[p],
                       parametric = F,
                       family = c("gaussian", "symmetric"),
                       method = c("loess", "model.frame"))
          
          # predict and compute NSE
          lwr.pred <- predict(lwr, y[ ,-ncol(y)])
          NS <- NSE(lwr.pred, y[ ,ncol(y)])
          if (NS > ns) {
            ns <- NS
            SelModel <- lwr
          }
        }
      }
      # Select the best four
    } else {
      cor.nw <- data.frame(idx = 1 : (ncol(x) - 1), cor = rep(0, ncol(x) - 1))
      for (g in 1:(ncol(x)-1)) {
        cor.nw[g, 2] <- ccf(x[ ,ncol(x)], x[ ,g], lag.max = 0, plot = F)$acf
      }
      cor.nw <- cor.nw[order(-cor.nw$cor),]
      best.lwr <- cor.nw$idx[1:4]
      
      # Set the dataframe for locally weighted regression
      lwr.tr <- x[ ,c(best.lwr, ncol(x))]
      lwr.te <- y[ ,c(best.lwr, ncol(y))]
      
      # Train the model
      for (g in 1:length(cands.lwr)) {
        for (p in 1:length(deg)) {
          lwr <- loess(out~., data = lwr.tr,  model = T,
                       span = cands.lwr[g], degree = deg[p],
                       parametric = F,
                       family = c("gaussian", "symmetric"),
                       method = c("loess", "model.frame"))
          
          # predict and compute NSE
          lwr.pred <- predict(lwr, lwr.te[ ,-ncol(lwr.te)])
          NS <- NSE(lwr.pred, lwr.te[ ,ncol(lwr.te)])
          if (NS > ns) {
            ns <- NS
            SelModel <- lwr
          }
        }
      }
    }
  }else {
    x <- as.data.frame(x)
    colnames(x)[ncol(x)] <- "out"
    dist <- c(1:23)
    ns <- -1000
    for (q in 1:length(dist)) {
      w.knn <- train.kknn(out~., data = x, kmax = 25, ks = NULL, distance = dist[q], kernel = c("rectangular", "triangular",  "optimal") )
      w.knn.pred <- predict(w.knn, y[ ,-ncol(y)])
      NS <- NSE(w.knn.pred, y[ ,ncol(y)])
      if (NS > ns) {
        ns <- NS
        SelModel <- w.knn
      }
    }
  }
  
  return(SelModel)
  
}

