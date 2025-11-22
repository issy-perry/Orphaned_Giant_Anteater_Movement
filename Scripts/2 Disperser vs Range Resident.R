# 2 Determining Who Is a Disperser vs Range-resident
# This script includes the following:
#     1. Check the DOF values for each wild-raised individual 
#     2. Check the DOF values for each orphaned individual 

# load packages
library(ctmm) # working with UDs

# import UDs
load("./RESULTS/AKDEs/UDs_orphan.rda") # orphaned 
load("./RESULTS/AKDEs/UDs_wild.rda") # wild-raised


# when DOF for area  on an AKDE is less than two, the individual is likely not range resident (https://groups.google.com/g/ctmm-user/c/R75QHDysv6s/m/rZX6BN4oAgAJ)
# looking at telemetry points and home ranges on a map is very informative on movement (individuals can just look like dispersers), but some individuals who look like they dispersed, may not have, 
# so quantifiying that with DOF evaluations is good practice



# wild-raised ----
summary(AKDE_wild[[1]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                            #312.6379 1607.5342 

summary(AKDE_wild[[2]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                          # 463.2915 2617.3976

summary(AKDE_wild[[3]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                          # 379.7297 2097.1895 

summary(AKDE_wild[[4]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                          # 312.6059 1577.6990

summary(AKDE_wild[[5]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                          # 386.2675 2152.8396 

summary(AKDE_wild[[6]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                        #   456.132  2622.395

summary(AKDE_wild[[7]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                        #   214.6436  913.7392 

summary(AKDE_wild[[8]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                        #   93.38042 283.43566 

summary(AKDE_wild[[9]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                        #   261.2825 1221.1188

summary(AKDE_wild[[10]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                          #   364.6994 1996.7240

summary(AKDE_wild[[11]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                        #   328.0711 1684.4358

summary(AKDE_wild[[12]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                          #   561.1717 3461.4349

summary(AKDE_wild[[13]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                          #   236.8982 1018.9096 

summary(AKDE_wild[[14]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                          #    163.3873  618.9248

summary(AKDE_wild[[15]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                        #     47.45445 109.25252

summary(AKDE_wild[[16]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                        #     302.0031 1504.4237

summary(AKDE_wild[[17]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                        #      231.4375 1001.3013 

summary(AKDE_wild[[18]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                        #      320.8466 1663.1603 

summary(AKDE_wild[[19]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                      #       429.9268 2451.1371

summary(AKDE_wild[[20]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                      #        301.4862 1573.6113

summary(AKDE_wild[[21]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                    #         309.3011 1545.0999 

summary(AKDE_wild[[22]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                    #         146.8317  523.0484  

summary(AKDE_wild[[23]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                      #         269.951  1217.914 

summary(AKDE_wild[[24]], units = FALSE) #disperser  $DOF area     bandwidth ********
                                              #         1.095371  1.061720 

summary(AKDE_wild[[25]], units = FALSE) #disperser  $DOF area     bandwidth ********
                                              #          1.153992  1.098126 

summary(AKDE_wild[[26]], units = FALSE) #disperser  $DOF area     bandwidth ********
                                            #           1.575748  1.490240 


# separate dispersers into lists
AKDE_wild_dis <- AKDE_wild[c(24,25,26)]
FITS_wild_dis <- FITS_wild[c(24,25,26)]

# separate range-residents into lists
AKDE_wild_RR <- AKDE_wild[-c(24,25,26)]
FITS_wild_RR <- FITS_wild[-c(24,25,26)]

# save outputs
save(AKDE_wild_dis, file = "./RESULTS/AKDEs/UDs_wild_dis.rda")
save(FITS_wild_dis, file =  "./RESULTS/Fits/Fits_wild_dis.rda")
save(AKDE_wild_RR, file = "./RESULTS/AKDEs/UDs_wild_RR.rda")
save(FITS_wild_RR, file = "./RESULTS/Fits/Fits_wild_RR.rda")




# orphaned ----
summary(AKDE_orphan[[1]], units = FALSE) #range-resident    $DOF area bandwidth 
                                                           #   108.0451  319.4879  

summary(AKDE_orphan[[2]], units = FALSE) #disperser    $DOF area bandwidth ***********
                                                      #   1.400545  1.416919

summary(AKDE_orphan[[3]], units = FALSE) #disperser    $DOF area bandwidth **********
                                                      #    1.294462  1.324629  

summary(AKDE_orphan[[4]], units = FALSE) #range-resident    $DOF area bandwidth 
                                                       #     236.0527  885.5536  

summary(AKDE_orphan[[5]], units = FALSE) #disperser    $DOF area bandwidth *********
                                                    #     1.970041  2.449087 

summary(AKDE_orphan[[6]], units = FALSE) #check again  $DOF area bandwidth *********
                                                  #     3.151928  3.556875 
#this ^^ is Cláudio, and he is definitely a disperser

summary(AKDE_orphan[[7]], units = FALSE) #range-resident    $DOF area bandwidth 
                                                         #    178.9701  580.8651  

summary(AKDE_orphan[[8]], units = FALSE) #range-resident     $DOF area bandwidth 
                                                        #     195.3897  644.0894 

summary(AKDE_orphan[[9]], units = FALSE) #range-resident     $DOF area bandwidth 
                                                      #       224.6922  815.8915 

summary(AKDE_orphan[[10]], units = FALSE) #range-resident      $DOF area bandwidth 
                                                        #       71.22162 184.44147 

summary(AKDE_orphan[[11]], units = FALSE) #disperser       $DOF area bandwidth **********
                                                      #      1.957936  1.962893  

summary(AKDE_orphan[[12]], units = FALSE) #range-resident $DOF area bandwidth 
                                                      #     74.75466 200.52628

summary(AKDE_orphan[[13]], units = FALSE) #range-resident    $DOF area bandwidth 
                                                          #     69.90854 187.78954 

summary(AKDE_orphan[[14]], units = FALSE) #disperser    $DOF area bandwidth *********
                                                     #      1.816522  2.105671 

summary(AKDE_orphan[[15]], units = FALSE) #check again     $DOF area bandwidth *********
                                                  #           4.874823  6.454469  
#this ^^ is Juju's first monitoring period, and she is definitely a disperser

summary(AKDE_orphan[[16]], units = FALSE) #range-resident    $DOF area bandwidth 
                                                       #       139.8215  421.4826  

summary(AKDE_orphan[[17]], units = FALSE) #range-resident $DOF area bandwidth 
                                                     #       16.70026  27.38263 

summary(AKDE_orphan[[18]], units = FALSE) #range-resident $DOF area bandwidth 
                                                    #      24.46184  44.74678

summary(AKDE_orphan[[19]], units = FALSE) #disperser  $DOF area bandwidth  *********
                                                      #       4.799933  6.006380  
#this ^^ is Nayeli and she has a clear dispersal period

summary(AKDE_orphan[[20]], units = FALSE) #range-resident $DOF area bandwidth 
                                                      #      13.18603  20.14422

summary(AKDE_orphan[[21]], units = FALSE) #disperser      $DOF area bandwidth **********
                                                    #        1.087802  1.099104

summary(AKDE_orphan[[22]], units = FALSE) #disperser        $DOF area bandwidth **********
                                                          #     1.344317  1.524567

summary(AKDE_orphan[[23]], units = FALSE)#disperser        $DOF area bandwidth **********
                                                         #     1.032718  1.050208

summary(AKDE_orphan[[24]], units = FALSE)#disperser        $DOF area bandwidth **********
                                                        #     1.576932  1.885042

summary(AKDE_orphan[[25]], units = FALSE)#range-resident        $DOF area bandwidth 
                                                              #     15.50106  24.88870

summary(AKDE_orphan[[26]], units = FALSE)#disperser        $DOF area bandwidth **********
                                                        #     1.242616  1.253588 

summary(AKDE_orphan[[27]], units = FALSE)#range-resident        $DOF area bandwidth 
                                                              #     28.43786  52.92036 

summary(AKDE_orphan[[28]], units = FALSE)#range-resident        $DOF area bandwidth 
                                                              #     40.88469  76.89346 


# based on the telemetry paths, there are some individuals who are either doing huge walkabouts or slight dispersal
# George, Heather, Mulan, and Peter (double check all)
# most of previous thoughts were true, but Juju is not a disperser as expected
# separate dispersers into lists
AKDE_orphan_dis <- AKDE_orphan[c(2,3,5,11,14,21,22,23,24,26)]
FITS_orphan_dis <- FITS_orphan[c(2,3,5,11,14,21,22,23,24,26)]

# separate range-residents into lists
AKDE_orphan_RR <- AKDE_orphan[-c(2,3,5,11,14,21,22,23,24,26)]
FITS_orphan_RR <- FITS_orphan[-c(2,3,5,11,14,21,22,23,24,26)]

# save outputs
save(AKDE_wild_dis, file = "./RESULTS/AKDEs/UDs_orphan_dis.rda")
save(FITS_wild_dis, file =  "./RESULTS/Fits/Fits_orphan_dis.rda")
save(AKDE_wild_RR, file = "./RESULTS/AKDEs/UDs_orphan_RR.rda")
save(FITS_wild_RR, file = "./RESULTS/Fits/Fits_orphan_RR.rda")




