# load packages
library(flextable)

# load data
Wild_info <- read.csv ("./Wild_info.csv")

# create tables
Wild_table <- flextable(Wild_info)

# increase line spacing so text is not so close to each other
Wild_table <- line_spacing(Wild_table, space = 1.5, part = "body")

# increase the widths of columns manually
Wild_table <- width(Wild_table, j = "ID", width = 2)
Wild_table <- width(Wild_table, j = "Sex", width = 1.5)
Wild_table <- width(Wild_table, j = "Age.Class", width = 1.5)
Wild_table <- width(Wild_table, j = "Monitoring.Period", width = 2)
Wild_table <- width(Wild_table, j = "Fixes", width = 1)

# save table
save_as_image(Wild_table, path = "./FIGURES/Revised/Supplementary_Table_2.png", 
              res = 600)

