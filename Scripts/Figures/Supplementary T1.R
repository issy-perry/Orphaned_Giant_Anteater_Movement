# load packages
library(flextable)

# load data
Orphan_info <- read.csv("./Orphan_info.csv")

# rewrite cells that loaded in wrong
Orphan_info[c(4), 3] <- "CRAS Barão de Cocais on 10/08/2023, Águas Vivas Retreat - Uberlândia on 02/21/2024" 
Orphan_info[c(5,9,12), 3] <- "TamanduASAS Support Center - Uberlândia" 

# create table
Orphan_table <- flextable(Orphan_info)

# increase line spacing so text is not so close to each other
Orphan_table <- line_spacing(Orphan_table, space = 1.5, part = "body")

# increase the widths of columns manually
Orphan_table <- width(Orphan_table, j = "ID", width = 0.5)
Orphan_table <- width(Orphan_table, j = "Sex", width = 0.5)
Orphan_table <- width(Orphan_table, j = "Rehabilitation.Center", width = 1.5)
Orphan_table <- width(Orphan_table, j = "Monitoring.Period", width = 1.5)
Orphan_table <- width(Orphan_table, j = "Fixes", width = 0.5)
Orphan_table <- width(Orphan_table, j = "Success", width = 1)
Orphan_table <- width(Orphan_table, j = "Context", width = 2.5)

# save table
save_as_image(Orphan_table, path = "./FIGURES/Revised/Supplementary_Table_1.png", 
              res = 600)

  