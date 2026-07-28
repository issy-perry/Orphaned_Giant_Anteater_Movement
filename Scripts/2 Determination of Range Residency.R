# 2 Determining Who Is a Disperser vs Range-resident
# This script includes the following:
#     1. Check the DOF values for each wild-raised individual 
#     2. Check the DOF values for each orphaned individual 
# when DOF for the area variable on an AKDE is less than two, the individual is likely not range resident (https://groups.google.com/g/ctmm-user/c/R75QHDysv6s/m/rZX6BN4oAgAJ)
# important to look at both the movement paths for evidence of dispersal events as well as the results of an AKDE.

# load packages
library(ctmm) # working with ctmm objects

# load results and data
# wild-raised
load("./RESULTS/AKDEs/UDs_wild.rda")
load("./RESULTS/Fits/Fits_wild.rda")
load("./RESULTS/Speed/Speed_wild.rda")
load("./DATA/Wild_raised_tel_data/Data_telemetry.rda")
# orphaned 
load("./RESULTS/Fits/Fits_orphan.rda") 
load("./RESULTS/AKDEs/UDs_orphan.rda")
load("./RESULTS/Speed/Speed_orphans.rda")
load("./DATA/Orphaned_tel_data/Data_telemetry.rda")

# wild-raised population ----
summary(AKDE_wild[["Alexander"]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                            #312.6379 1607.5342 

summary(AKDE_wild[["Annie"]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                          # 463.2915 2617.3976

summary(AKDE_wild[["Anthony"]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                          # 379.7297 2097.1895 

summary(AKDE_wild[["Beto"]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                          # 312.6059 1577.6990

summary(AKDE_wild[["Bumpus"]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                          # 386.2675 2152.8396 

summary(AKDE_wild[["Cate"]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                        #   456.132  2622.395

summary(AKDE_wild[["Christoffer"]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                        #   214.6436  913.7392 

summary(AKDE_wild[["Elaine"]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                        #   93.38042 283.43566 

summary(AKDE_wild[["Hannah"]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                        #   261.2825 1221.1188

summary(AKDE_wild[["Jackson"]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                          #   364.6994 1996.7240

summary(AKDE_wild[["Jane"]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                        #   328.0711 1684.4358

summary(AKDE_wild[["Kyle"]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                          #   561.1717 3461.4349

summary(AKDE_wild[["Larry"]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                          #   236.8982 1018.9096 

summary(AKDE_wild[["Little_Rick"]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                          #    163.3873  618.9248

summary(AKDE_wild[["Luigi"]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                        #     47.45445 109.25252

summary(AKDE_wild[["Makao"]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                        #     302.0031 1504.4237

summary(AKDE_wild[["Margaret"]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                        #      231.4375 1001.3013 

summary(AKDE_wild[["Maria"]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                        #      320.8466 1663.1603 

summary(AKDE_wild[["Puji"]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                      #       429.9268 2451.1371

summary(AKDE_wild[["Reid"]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                      #        301.4862 1573.6113

summary(AKDE_wild[["Rodolfo"]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                    #         309.3011 1545.0999 

summary(AKDE_wild[["Sheron"]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                    #         146.8317  523.0484  

summary(AKDE_wild[["Thomas"]], units = FALSE) #range-resident  $DOF area     bandwidth 
                                                      #         269.951  1217.914 

summary(AKDE_wild[["Delphine"]], units = FALSE) #disperser  $DOF area     bandwidth ********
                                              #         1.095371  1.061720 

summary(AKDE_wild[["Gala"]], units = FALSE) #disperser  $DOF area     bandwidth ********
                                              #          1.153992  1.098126 

summary(AKDE_wild[["Segre"]], units = FALSE) #disperser  $DOF area     bandwidth ********
                                            #           1.575748  1.490240 

# separate range-residents into lists
AKDE_wild <- AKDE_wild[c("Alexander", "Annie", "Anthony", "Beto", "Bumpus", "Cate", "Christoffer", 
                         "Elaine", "Hannah", "Jackson", "Jane", "Kyle", "Larry", "Little_Rick", 
                         "Luigi", "Makao", "Margaret", "Maria", "Puji", "Reid", "Rodolfo", "Sheron", 
                         "Thomas")] 
FITS_wild <- FITS_wild[c("Alexander", "Annie", "Anthony", "Beto", "Bumpus", "Cate", "Christoffer", 
                         "Elaine", "Hannah", "Jackson", "Jane", "Kyle", "Larry", "Little_Rick", 
                         "Luigi", "Makao", "Margaret", "Maria", "Puji", "Reid", "Rodolfo", "Sheron", 
                         "Thomas")] 
SPEED_wild <- SPEED_wild[c("Alexander", "Annie", "Anthony", "Beto", "Bumpus", "Cate", "Christoffer", 
                           "Elaine", "Hannah", "Jackson", "Jane", "Kyle", "Larry", "Little_Rick", 
                           "Luigi", "Makao", "Margaret", "Maria", "Puji", "Reid", "Rodolfo", "Sheron", 
                           "Thomas")]
DATA_wild <- DATA_wild[c("Alexander", "Annie", "Anthony", "Beto", "Bumpus", "Cate", "Christoffer", 
                         "Elaine", "Hannah", "Jackson", "Jane", "Kyle", "Larry", "Little_Rick", 
                         "Luigi", "Makao", "Margaret", "Maria", "Puji", "Reid", "Rodolfo", "Sheron", 
                         "Thomas")]
# save outputs
save(AKDE_wild, file = "./RESULTS/AKDEs/UDs_wild_RR.rda")
save(FITS_wild, file = "./RESULTS/Fits/Fits_wild_RR.rda")
save(SPEED_wild, file = "./RESULTS/Speed/Speed_wild_RR.rda")
save(DATA_wild, file = "./DATA/Wild_raised_tel_data/Data_telemetry_RR.rda")

# orphaned population ----
summary(AKDE_orphan[["Arya"]], units = FALSE) #range-resident    $DOF area bandwidth 
                                                           #   108.0451  319.4879  

summary(AKDE_orphan[["Bahia"]], units = FALSE) #disperser    $DOF area bandwidth ***********
                                                      #   1.400545  1.416919

summary(AKDE_orphan[["Beezie"]], units = FALSE) #disperser    $DOF area bandwidth **********
                                                      #    1.294462  1.324629  

summary(AKDE_orphan[["Bella"]], units = FALSE) #range-resident    $DOF area bandwidth 
                                                       #     236.0527  885.5536  

summary(AKDE_orphan[["Capitu"]], units = FALSE) #disperser    $DOF area bandwidth *********
                                                    #     1.970041  2.449087 

summary(AKDE_orphan[["Cláudio"]], units = FALSE) #check again  $DOF area bandwidth *********
                                                  #     3.151928  3.556875 
#this ^^ is definitely a disperser

summary(AKDE_orphan[["Colete"]], units = FALSE) #range-resident    $DOF area bandwidth 
                                                         #    178.9701  580.8651  

summary(AKDE_orphan[["Dom"]], units = FALSE) #range-resident     $DOF area bandwidth 
                                                        #     195.3897  644.0894 

summary(AKDE_orphan[["Dumbo_1"]], units = FALSE) #range-resident     $DOF area bandwidth 
                                                      #       224.6922  815.8915 

summary(AKDE_orphan[["Dumbo_2"]], units = FALSE) #range-resident      $DOF area bandwidth 
                                                        #       71.22162 184.44147 

summary(AKDE_orphan[["Erick"]], units = FALSE) #disperser       $DOF area bandwidth **********
                                                      #      1.957936  1.962893  

summary(AKDE_orphan[["George"]], units = FALSE) #range-resident $DOF area bandwidth 
                                                      #     74.75466 200.52628

summary(AKDE_orphan[["Heather"]], units = FALSE) #range-resident    $DOF area bandwidth 
                                                          #     69.90854 187.78954 

summary(AKDE_orphan[["Jacobina"]], units = FALSE) #disperser    $DOF area bandwidth *********
                                                     #      1.816522  2.105671 

summary(AKDE_orphan[["Juju_1"]], units = FALSE) #check again     $DOF area bandwidth *********
                                                  #           4.874823  6.454469  
# definitely a disperser

summary(AKDE_orphan[["Juju_2"]], units = FALSE) #range-resident    $DOF area bandwidth 
                                                       #       139.8215  421.4826  

summary(AKDE_orphan[["Mulan"]], units = FALSE) #range-resident $DOF area bandwidth 
                                                     #       16.70026  27.38263 

summary(AKDE_orphan[["Nancy"]], units = FALSE) #range-resident $DOF area bandwidth 
                                                    #      24.46184  44.74678

summary(AKDE_orphan[["Nayeli"]], units = FALSE) #disperser  $DOF area bandwidth  *********
                                                      #       4.799933  6.006380  
#this ^^ is Nayeli and she has a clear dispersal period

summary(AKDE_orphan[["Peter"]], units = FALSE) #range-resident $DOF area bandwidth 
                                                      #      13.18603  20.14422

summary(AKDE_orphan[["Renee_1"]], units = FALSE) #disperser      $DOF area bandwidth **********
                                                    #        1.087802  1.099104

summary(AKDE_orphan[["Renee_2"]], units = FALSE) #disperser        $DOF area bandwidth **********
                                                          #     1.344317  1.524567

summary(AKDE_orphan[["Renee_3"]], units = FALSE)#disperser        $DOF area bandwidth **********
                                                         #     1.032718  1.050208

summary(AKDE_orphan[["Renee_4"]], units = FALSE)#disperser        $DOF area bandwidth **********
                                                        #     1.576932  1.885042

summary(AKDE_orphan[["Rita"]], units = FALSE)#range-resident        $DOF area bandwidth 
                                                              #     15.50106  24.88870

summary(AKDE_orphan[["Tim_1"]], units = FALSE)#disperser        $DOF area bandwidth **********
                                                        #     1.242616  1.253588 

summary(AKDE_orphan[["Tim_2"]], units = FALSE)#range-resident        $DOF area bandwidth 
                                                              #     28.43786  52.92036 

summary(AKDE_orphan[["Tim_3"]], units = FALSE)#range-resident        $DOF area bandwidth 
                                                              #     40.88469  76.89346 


# remove dispersers from dataset
AKDE_orphan <- AKDE_orphan[c("Arya", "Bella", "Colete", "Dom", "Dumbo_1", "Dumbo_2", "George", 
                             "Heather", "Juju_2", "Mulan", "Nancy", "Peter", "Rita", "Tim_2", 
                             "Tim_3")]
FITS_orphan <- FITS_orphan[c("Arya", "Bella", "Colete", "Dom", "Dumbo_1", "Dumbo_2", "George", 
                             "Heather", "Juju_2", "Mulan", "Nancy", "Peter", "Rita", "Tim_2", 
                             "Tim_3")] 
SPEED_orphan <- SPEED_orphan[c("Arya", "Bella", "Colete", "Dom", "Dumbo_1", "Dumbo_2", "George", 
                               "Heather", "Juju_2", "Mulan", "Nancy", "Peter", "Rita", "Tim_2", 
                               "Tim_3")]
DATA_orphan <- DATA_orphan[c("Arya", "Bella", "Colete", "Dom", "Dumbo_1", "Dumbo_2", "George", 
                             "Heather", "Juju_2", "Mulan", "Nancy", "Peter", "Rita", "Tim_2", 
                             "Tim_3")]
# save outputs
save(AKDE_orphan, file = "./RESULTS/AKDEs/UDs_orphan_RR.rda")
save(FITS_orphan, file = "./RESULTS/Fits/Fits_orphan_RR.rda")
save(SPEED_orphan, file = "./RESULTS/Speed/Speed_orphans_RR.rda")
save(DATA_orphan, file = "./DATA/Orphaned_tel_data/Data_telemetry_RR.rda")
