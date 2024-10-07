# this is the script for cleaning the data for the modeling, before running run master function list
yaml_data <- yaml.load_file("./Code/config_site.yml")
clades <- yaml_data$clades # the clades for the week
file_name <- paste(yaml_data$data, yaml_data$target_date)
file_name <- paste0(file_name, ".tsv.gz")
target_date <- as.Date(yaml_data$target_date)
data <- read_tsv(file_name)
data_used <- trim_clades(data = data, clades = clades)
data_used <- filter(data_used, location %in% c(state.name, "Puerto Rico", "Washington DC"))# removing 
# locations we don't want to model