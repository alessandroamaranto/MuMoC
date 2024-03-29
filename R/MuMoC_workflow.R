# MuMoC workflow: a multimodel combination procedure for improving GW forecasts
# in the High-Plains aquifer (USA). MuMoC is composed by 3 main components as
# follows: 
# 1 - ANN routine: a multilayer perceptron is run for each well, and forecasts 
#                  GW level up to four months ahead. Forecasts having accuracy
#                  lower than 0.6 in NSE are classified as 'bad'.
# 2-  WS routine:  for each of the 'bad' well, the geographical neighborhood is
#                  inspected. The wells maximizing the selection criteria (pareto
#                  -optimal) are selected.
# 3- MC routine:   Forecasts from the neighborhood are combined with IBL, and per
#                  formances tested

# libraries
source("MuMoC_fun.R")
library(rgeos)
library(sp)
library(geosphere)
library(caret)
library(foreach)
library(doParallel)
library(rPref)
library(ggplot2)


## Import data
load('Data_in/MuMoC_in.RData')
load('Data_in/candidate.RData')

# load column index and river-excluded text files
c <- scan("Data_in/column_index.txt", numeric())
exc <- scan("Data_in/river_out.txt", logical()) 

# import flow data and metadata from USGS 
flow_m <- load_flow_s('Data_in/hp.shp', "00060", "1980-01-01", "2017-08-01")
flow <- load_flow(flow_m, "00060", "1980-01-01", "2017-11-30")

# exclude non-selected rivers
flow_m <- exclude_meta(flow_m, flow)
riv_idx <- min_distance(s[ ,c(5, 4)], coordinates(flow_m))

# initialize
flow_c <- h <- list()
N <- length(X); Nh <- 4
Vs <- d_eps <- eps_v <- matrix(0, N, Nh)

# ANN routine (warning: combined with IVS might take several days)

for (i in 1:N) {
  
  Nvar <- ncol(X[[i]])
  flow_c[[i]] <- river_coupling(X[[i]], flow[[riv_idx[i, 1]]], exc[i])
  Xc <- normalize(X[[i]][, -1], split[[i]][, 2])
  Xv <- Xc[[2]]
  Xc <- Xc[[1]]
  
  XRc <- add_river(flow_c[[i]][,-1], split[[i]][, 2], exc[i])
  if(!is.character(XRc)){
    XRv <- XRc[[2]]
    XRc <- XRc[[1]]
  }
  
  Xc <- attach_date(X[[i]], Xc, split[[i]][, 2], 'calibration')
  Xv <- attach_date(X[[i]], Xv, split[[i]][, 2], 'validation')
  Fl <- createFolds(Xc[, Nvar], 10, list = T)
  
  for (j in 1:Nh) {
    
    xc <- select_columns(j, Xc, c)
    xv <- select_columns(j, Xv, c)
    
    Ir <- input_variable_selection(cand, xc, xv, e, Fl, XRc, XRv, j, 'run')
    Vs[i, j] <- Ir[[1]]
    d_eps[i, j] <- Ir[[2]]
    
    eps_v[i, j] <- ann_err(cand, xc, xv, e, Fl, XRc, XRv, j, Vs[i, j])     
  }
}

# Results, section 1, figures and tables
eps_state <- table_3(eps_v, s, 'Data_in/USA.shp')

# MuMoC starts: identify 'b' models, get the opt neigh
b <- apply(eps_v, 2, function(x) which(x < 0.6))
Pf <- exhaustive_pareto_search(b, X, s, eps_v)

# Models combination-------------------------------

# Prepare data for models combination
MuMoC_in <- prepare_MuMoC_input(b, eps_v, Pf, Vs, riv_idx, X, split, exc, flow, c)
Xc <- MuMoC_in[[1]]; Xv <- MuMoC_in[[2]]
eps_MuMoC <- list()

# Run the combination
for (k in 1:Nh) {
  
  Xck <- Xc[[k]]
  Xvk <- Xv[[k]]
  
  eps0 <- return_IBL_error(Xck, Xvk)  
  eps0[ ,6] <- eps_v[b[[k]], k]
  eps_MuMoC[[k]] <- eps0 
}

# Plot the results------------------

# Allocate memory for plot database
ssel <- list()

# Iterate along lead times
for (k in 1:4) {
  
  # Initialize pareto list
  par.ls <- pareto.lt[[k]]
  sel <- matrix(nrow = length(par.ls), ncol = 4)
  
  # Iterate along wells
  for (i in 1:length(par.ls)) {
    
    # Allocate memory for distance vector
    par <- par.ls[[i]]
    d <- vector(length = nrow(par))
    
    # Compute distance
    for (j in 1:nrow(par)) {
      d[j] <- sqrt( (par$NSE[j]-1)^2 + (par$Corr[j]-1)^2)
    }
    
    # Select the right pareto point
    sel[i, ] <- as.numeric(c(par[which.min(d), c(1, 2)], res.tot[[k]][i ,3], res.tot[[k]][i ,3] - res.tot[[k]][i ,6]))
  }
  
  
  if(k == 3){
    ssel[[k]] <- as.data.frame(sel[-28, ] )
  }else if(k == 4){
    ssel[[k]] <- as.data.frame(sel[-38, ] )
  }else{
    ssel[[k]] <- as.data.frame(sel)
  }
}

new.pareto.lt <- pareto.lt
# Perform the same operation with the updated values on the correct training and testing
for (k in 1:n.steps) {
  
  # Select the bad models and the code for the neighbours
  bad.models <- which(res.nse[ ,k] < 0.6)
  
  train.well <- list()
  test.well <- list()
  
  alternatives <- matrix(nrow = length(bad.models), ncol = 4)
  for (i in 1:length(bad.models)) {
    
    print(paste0(k, "-", i))
    neigh.code <- pareto.lt[[k]][[i]]$Wid
    
    if (sum(neigh.code) != 0){
      
      mod.train <- ibl.train[[k]][[i]]
      mod.test <- ibl.test[[k]][[i]]
      
      # Select neighbouring wells from the list, as well as their ivs results
      neigh.list <- ivs.eff[neigh.code]
      n.iv <- ns.ivs[neigh.code]
      
      # Select dates for split
      split.neigh <- data.frame(Date = as.data.frame(split.eff[  bad.models[i]  ])[, 1])
      split<- split.eff[[bad.models[i]]][ ,2]
      
      # Select well data
      well <- ivs.eff[[bad.models[i]]]
      
      # Divide train and test set
      well.fake <- list(well)
      
      # Return them
      train.neigh <- mapply(NeighTrainTest, neigh.list, split.neigh, well.fake, 1, SIMPLIFY = F) 
      test.neigh <- mapply(NeighTrainTest, neigh.list, split.neigh, well.fake, 2, SIMPLIFY = F) 
      
      # Select only WL column
      train <- lapply(train.neigh, function(x) y <- as.data.frame(x[ ,43]) )
      train <- do.call(cbind, train)
      
      test <- lapply(test.neigh, function(x) y <- as.data.frame(x[ ,43]) )
      test <- do.call(cbind, test)
      
      # Normalize
      max.par <- apply(train, 2, max)
      min.par <- apply(train, 2, min)
      
      
      for (q in 1:ncol(train)) {
        train[ ,q] <- (train[ ,q] - min.par[q])/(max.par[q] - min.par[q])
        test[ ,q] <- (test[ ,q] - min.par[q])/(max.par[q] - min.par[q])
        new.pareto.lt[[k]][[i]][q, 1] <- NSE(test[ ,q], mod.test[ ,q] )
        fakecorr.in <- c(train[ ,q], test[ ,q])
        fakecorr.out <- c(mod.train[ ,ncol(mod.train)], mod.test[ ,ncol(mod.test)])
        new.pareto.lt[[k]][[i]][q, 2] <- ccf(fakecorr.in, fakecorr.out, plot = F, lag.max = 0)$acf
      }
    }
  }
}


ssel.2 <- list()

# Iterate along lead times
for (k in 1:4) {
  
  # Initialize pareto list
  par.ls <- new.pareto.lt[[k]]
  sel <- matrix(nrow = length(par.ls), ncol = 4)
  
  # Iterate along wells
  for (i in 1:length(par.ls)) {
    
    # Allocate memory for distance vector
    par <- par.ls[[i]]
    d <- vector(length = nrow(par))
    
    # Compute distance
    for (j in 1:nrow(par)) {
      d[j] <- sqrt( (par$NSE[j]-1)^2 + (par$Corr[j]-1)^2)
    }
    
    # Select the right pareto point
    sel[i, ] <- as.numeric(c(par[which.min(d), c(1, 2)], res.tot[[k]][i ,3], res.tot[[k]][i ,3] - res.tot[[k]][i ,6]))
  }
  
  
  if(k == 3){
    ssel.2[[k]] <- as.data.frame(sel[-28, ] )
  }else if(k == 4){
    ssel.2[[k]] <- as.data.frame(sel[-38, ] )
  }else{
    ssel.2[[k]] <- as.data.frame(sel)
  }
}

# Plot results: Objective 3--------------------------------

# Objective 3.2: Kansas and Nebraska
hp.states <- c("Kansas", "Nebraska")
us.used <- us.used[which(us.used@data$STATE_NAME %in% hp.states), ]

# Allocate the memory for the indicator
new.nse.obj3 <- matrix(nrow = length(hp.states), ncol = 4)
old.nse.obj3 <- matrix(nrow = length(hp.states), ncol = 4)
shift.perc <- matrix(nrow = 2, ncol = 4)
not.available <- matrix(nrow = 2, ncol =4)
# Iterate along the states
for (i in 1:length(hp.states)) {
  
  # Extract ith State
  state <- us.used[which(us.used@data$STATE_NAME %in% hp.states[i]), ]
  
  for (k in 1:n.steps) {
    
    # Extract bad models
    res <- res.tot[[k]][ ,c(3, 6)]
    bad.models.location <- new.meta[which(res.nse[ ,k] < 0.6), ]
    
    # Remove weird 38th well
    if (k == 3 ) {
      res = res[-28, ]
      bad.models.location <- bad.models.location[-28, ]
    }else if (k == 4) {
      bad.models.location <- bad.models.location[-38, ]
      res = res[-38, ]
    }
    
    # Take data from the state
    points <- which(!is.na (over(bad.models.location, state,  fn = NULL))[ ,1]) 
    res <- res[points, ]
    bad.models.location <- bad.models.location[points, ]
    
    # Remove missing values
    remove <- which(is.na(res[, 1]))
    
    if (length(remove) > 0){
      res <- res[-remove, ]
      bad.models.location <- bad.models.location[-remove, ]
    }
    not.available[i, k] <- length(remove)
    new.nse.obj3[i ,k] <- mean(res[ ,1])
    old.nse.obj3[i ,k] <- mean(res[ ,2])
    shift.perc[i, k] <- length(which(res[, 1] > 0.6))/length(res[, 1])
  }
  
}

# Manipulate a little to get the data in the desired form
kansas.delta.nse <- data.frame(LT = c(1:4), knn = new.nse.obj3[1, ], ann = old.nse.obj3[1, ])
nebraska.delta.nse <- data.frame(LT = c(1:4), knn = new.nse.obj3[2, ], ann = old.nse.obj3[2, ])

# Reshape
kansas.delta.nse <- melt(kansas.delta.nse, "LT")
nebraska.delta.nse <- melt(nebraska.delta.nse, "LT")

kansas.shift <- data.frame(LT = c(1:4), shift = shift.perc[1, ])
nebraska.shift <- data.frame(LT = c(1:4), shift = shift.perc[2, ])

kansas.shift <- melt(kansas.shift, "LT")
nebraska.shift <- melt(nebraska.shift, "LT")
# Plot

p.ks <- ggplot(kansas.delta.nse, aes(LT, value, group = variable, col = variable))+
  # geom_line()+
  geom_point(size = 2)+
  ylab(expression(mu ~ "(NSE)"))+
  xlab(expression("Lead Time (months)"))+
  ggtitle("Kansas")+
  scale_y_continuous(limits = c(0.38, 0.58))+
  theme_bw()+
  #scale_colour_manual(values = c("cyan", "purple" , "cyan", "magenta", "magenta")) +  
  #scale_linetype_manual(values=c("dashed", "longdash", "solid", "solid","dashed")) + 
  #scale_size_manual(values = c(1.2, 1.2, 1.2, 1.2, 1.2))+
  theme(text=element_text(family="Times New Roman", size=12),
        legend.key.size = unit(2,"line"), legend.text=element_text(size=12))+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        plot.margin = unit(c(1, 2, 0, 0 ), "pt" ) )


p.ne <- ggplot(nebraska.delta.nse, aes(LT, value, group = variable, col = variable))+
  # geom_line()+
  geom_point(size = 2)+
  ylab("")+
  xlab(expression("Lead Time (months)"))+
  ggtitle("Nebraska")+
  scale_y_continuous(limits = c(0.38, 0.58))+
  theme_bw()+
  #scale_colour_manual(values = c("cyan", "purple" , "cyan", "magenta", "magenta")) +  
  #scale_linetype_manual(values=c("dashed", "longdash", "solid", "solid","dashed")) + 
  #scale_size_manual(values = c(1.2, 1.2, 1.2, 1.2, 1.2))+
  theme(text=element_text(family="Times New Roman", size=12),
        legend.key.size = unit(2,"line"), legend.text=element_text(size=12))+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        plot.margin = unit(c(1, 0, 0, 2 ), "pt" )     )

sh.kans <- ggplot(kansas.shift, aes(LT, value, group = variable, col = variable))+
  # geom_line()+
  geom_point(size = 2)+
  ylab(expression("% of class switch"))+
  xlab(expression(""))+
  scale_y_reverse(limits = c(0.6, 0.0))+
  theme_bw()+
  #scale_colour_manual(values = c("cyan", "purple" , "cyan", "magenta", "magenta")) +  
  #scale_linetype_manual(values=c("dashed", "longdash", "solid", "solid","dashed")) + 
  #scale_size_manual(values = c(1.2, 1.2, 1.2, 1.2, 1.2))+
  theme(text=element_text(family="Times New Roman", size=12),
        legend.key.size = unit(2,"line"), legend.text=element_text(size=12))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        plot.margin = unit(c(0, 2, 0, 0 ), "pt" ))

sh.ne <- ggplot(nebraska.shift, aes(LT, value, group = variable, col = variable))+
  # geom_line()+
  geom_point(size = 2)+
  ylab("")+
  xlab(expression(""))+
  scale_y_reverse(limits = c(0.6, 0.0))+
  theme_bw()+
  #scale_colour_manual(values = c("cyan", "purple" , "cyan", "magenta", "magenta")) +  
  #scale_linetype_manual(values=c("dashed", "longdash", "solid", "solid","dashed")) + 
  #scale_size_manual(values = c(1.2, 1.2, 1.2, 1.2, 1.2))+
  theme(text=element_text(family="Times New Roman", size=12),
        legend.key.size = unit(2,"line"), legend.text=element_text(size=12))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        plot.margin = unit(c(0, 0, 0, 1 ), "pt" ))

p.legend <- ggplot(nebraska.delta.nse, aes(LT, value, group = variable, col = variable))+
  # geom_line()+
  geom_point(size = 2)+
  ylab(expression(mu ~ "(NSE)"))+
  xlab(expression("Lead Time (months)"))+
  scale_y_continuous(limits = c(0.38, 0.57))+
  theme_bw()+
  scale_color_manual(values = cols,  labels = c(expression(italic(MMC)),
                                                expression(italic(Ann))))+
  labs(col = "Model:")

library(cowplot)

prow <- plot_grid( p.ks + theme(legend.position="none"),
                   p.ne + theme(legend.position="none"),
                   sh.kans + theme(legend.position="none"),
                   sh.ne + theme(legend.position="none"),
                   align = 'hv',
                   label_size = 9,
                   hjust = -1,
                   nrow = 2
)


prow

legend_b <- get_legend(p.legend+ theme(legend.position = "bottom"))
p5 <- plot_grid( prow, legend_b, ncol = 1, rel_heights  = c(0.7, 0.03))
p5

ggsave("R32.tiff", plot = last_plot(), device = "tiff", path = NULL, scale = 1, dpi = 600, limitsize = T)
# Objective 3.3------------------------------------
# Plot element one
x <- ssel.2[[1]]

p.plot.1 <- ggplot(data = x, aes(V1, V2, size = V3, col = V4))+
  geom_point()+
  xlab("Neighbour Well NSE")+
  ylab("Cross Correlation")+
  ggtitle("Lead Time = 1 Month"   )+
  scale_color_gradient2(low = "white", mid = "dodgerblue", high = "dodgerblue4",
                        midpoint = 0.0,
                        #breaks = c(quantile(aaa$V4, c(0.01, 0.10, 0.50, 0.75, 0.90))),
                        #breaks = round(seq (min(x$V4), max(x$V4),  (max(x$V4) - min(x$V4)) / 4  ),2) ,
                        breaks = c(-0.05, 0,  0.05, 0.1) ,
                        limits = c(min(x$V4), max(x$V4 ) ) )+
  #scale_size_continuous(breaks = round(as.numeric(quantile(x$V3, c(0.01, 0.25, 0.5, 0.75, 0.9)) ),2 ), range = c(3, 12))+
  scale_size_continuous(breaks = c(0.35, 0.5, 0.6, 0.65), range = c(3, 12))+
  theme_bw()+
  theme(text = element_text(family = "Times New Roman",  size = 10))+
  theme(axis.text=element_text(family = "Times New Roman", size = 10),
        axis.title=element_text(family = "Times New Roman", size=12, face="bold"), 
        legend.title=element_text(family = "Times New Roman", size=12, face = "bold"), 
        legend.text = element_text(margin = margin(r = 28, unit = "pt")), 
        plot.title = element_text(family = "Times New Roman", size = 14, face = "bold", hjust = 0.5), 
        plot.margin = margin(0.1 ,0.5, 0.5, 0, unit = "cm") )+
  labs(col = expression(paste(Delta, "NSE")), size = "NSE")
#p.plot.1


x <- ssel.2[[2]]

p.plot.2 <- ggplot(data = x, aes(V1, V2, size = V3, col = V4))+
  geom_point()+
  xlab("Neighbour Well NSE")+
  ylab("Cross Correlation")+
  ggtitle("Lead Time = 2 Months"   )+
  scale_color_gradient2(low = "aliceblue", mid = "dodgerblue", high = "dodgerblue4",
                        midpoint = 0.0,
                        #breaks = c(quantile(aaa$V4, c(0.01, 0.10, 0.50, 0.75, 0.90))),
                        #breaks = round(seq (min(x$V4), max(x$V4),  (max(x$V4) - min(x$V4)) / 4  ),2) ,
                        breaks = c(0, 0.1, 0.2) ,
                        limits = c(min(x$V4), max(x$V4 ) ) )+
  #scale_size_continuous(breaks = round(as.numeric(quantile(x$V3, c(0.01, 0.25, 0.5, 0.75, 0.9)) ),2 ), range = c(3, 12))+
  scale_size_continuous(breaks = c(0.30, 0.50, 0.60, 0.70), range = c(3, 12))+
  theme_bw()+
  theme(text = element_text(family = "Times New Roman",  size = 10))+
  theme(axis.text=element_text(family = "Times New Roman", size = 10),
        axis.title=element_text(family = "Times New Roman", size=12, face="bold"), 
        legend.title=element_text(family = "Times New Roman", size=12, face = "bold"), 
        legend.text = element_text(margin = margin(r = 28, unit = "pt")), 
        plot.title = element_text(family = "Times New Roman", size = 14, face = "bold", hjust = 0.5), 
        plot.margin = margin(0.1 ,0.1, 0.5, 0.5, unit = "cm") )+
  labs(col = expression(paste(Delta, "NSE")), size = "NSE")+
  guides(color = guide_colorbar(order = 1),
         size = guide_legend(order = 0)
         
  )

p.plot.2


x <- ssel.2[[3]]

p.plot.3 <- ggplot(data = x, aes(V1, V2, size = V3, col = V4))+
  geom_point()+
  #scale_x_continuous(limits = c(0.5, 1))+
  #scale_y_continuous(limits = c(0.5, 1))+
  xlab("Neighbour Well NSE")+
  ylab("Cross Correlation")+
  ggtitle("Lead Time = 3 Months"   )+
  scale_color_gradient2(low = "aliceblue", mid = "dodgerblue", high = "dodgerblue4",
                        midpoint = 0.0,
                        #breaks = c(quantile(aaa$V4, c(0.01, 0.10, 0.50, 0.75, 0.90))),
                        #breaks = round(seq (min(x$V4), max(x$V4),  (max(x$V4) - min(x$V4)) / 4  ),2) ,
                        breaks = c(-0.2, -0.1, 0, 0.1, 0.2) ,
                        limits = c(min(x$V4), max(x$V4 ) ) )+
  #scale_size_continuous(breaks = round(as.numeric(quantile(x$V3, c(0.01, 0.25, 0.5, 0.75, 0.9)) ),2 ), range = c(3, 12))+
  scale_size_continuous(breaks = c(0.20, 0.50, 0.60, 0.70), range = c(3, 12))+
  theme_bw()+
  theme(text = element_text(family = "Times New Roman",  size = 10))+
  theme(axis.text=element_text(family = "Times New Roman", size = 10),
        axis.title=element_text(family = "Times New Roman", size=12, face="bold"), 
        legend.title=element_text(family = "Times New Roman", size=12, face = "bold"), 
        legend.text = element_text(margin = margin(r = 28, unit = "pt")), 
        plot.title = element_text(family = "Times New Roman", size = 14, face = "bold", hjust = 0.5), 
        plot.margin = margin(0.5 ,0.5, 0.5, 0.0, unit = "cm") )+
  labs(col = expression(paste(Delta, "NSE")), size = "NSE")+
  guides(color = guide_colorbar(order = 1),
         size = guide_legend(order = 0)
         
  )

x <- ssel.2[[4]]

p.plot.4 <- ggplot(data = x, aes(V1, V2, size = V3, col = V4))+
  geom_point()+
  #scale_x_continuous(limits = c(0.5, 1))+
  #scale_y_continuous(limits = c(0.5, 1))+
  xlab("Neighbour Well NSE")+
  ylab("Cross Correlation")+
  ggtitle("Lead Time = 4 Months"   )+
  scale_color_gradient2(low = "aliceblue", mid = "dodgerblue", high = "dodgerblue4",
                        midpoint = 0.0,
                        #breaks = c(quantile(aaa$V4, c(0.01, 0.10, 0.50, 0.75, 0.90))),
                        #breaks = round(seq (min(x$V4), max(x$V4),  (max(x$V4) - min(x$V4)) / 4  ),2) ,
                        breaks = c(-0.2, -0.1, 0, 0.1, 0.2, 0.3) ,
                        limits = c(min(x$V4), max(x$V4 ) ) )+
  #scale_size_continuous(breaks = round(as.numeric(quantile(x$V3, c(0.01, 0.25, 0.5, 0.75, 0.9)) ),2 ), range = c(3, 12))+
  scale_size_continuous(breaks = c(0.25, 0.50, 0.60, 0.70), range = c(3, 12))+
  theme_bw()+
  theme(text = element_text(family = "Times New Roman",  size = 10))+
  theme(axis.text=element_text(family = "Times New Roman", size = 10),
        axis.title=element_text(family = "Times New Roman", size=12, face="bold"), 
        legend.title=element_text(family = "Times New Roman", size=12, face = "bold"), 
        legend.text = element_text(margin = margin(r = 28, unit = "pt")), 
        plot.title = element_text(family = "Times New Roman", size = 14, face = "bold", hjust = 0.5), 
        plot.margin = margin(0.5 ,0.1, 0.5, 0.5, unit = "cm") )+
  labs(col = expression(paste(Delta, "NSE")), size = "NSE")+
  guides(color = guide_colorbar(order = 1),
         size = guide_legend(order = 0)
         
  )
p.plot.4

library(cowplot)
lt.1.4 <- plot_grid(p.plot.1,
                    p.plot.2,
                    p.plot.3,
                    p.plot.4,
                    align = "vh",
                    labels = c("A", "B", "C", "D"),
                    label_size = 9,
                    hjust = -1, 
                    nrow = 2)

ggsave("Obj3.tiff", plot = last_plot(), device = "tiff", path = NULL, scale = 1, dpi = 600, limitsize = T)

x <- do.call(rbind, ssel.2)

