How To Directly Access OCO-2 Data from an S3 Bucket using R
================

### Author: Jon Hobbs, Alexis Hunzinger

### Date Authored: 2025-08-12

## Overview

This notebook provides a quick demonstration on accessing and
summarizing products from the Orbiting Carbon Observatory-2 (OCO-2)
hosted via an Amazon S3 bucket. It demonstrates how to access an S3
bucket with the `aws.s3` package, how to read data with the `hdf5r`
package, and how to quickly compute and plot small-area aggregate
estimates of atmospheric carbon dioxide using data wrangling tools from
the `tidyverse` and visualization with `ggplot2`.

The [OCO-2](https://doi.org/10.5067/70K2B2W8MNGY) and
[OCO-3](https://doi.org/10.5067/8U0VGVQC7HZG) missions provide estimates
of atmospheric carbon dioxide (XCO2) and solar-induced fluorescence
(SIF) at fine spatial resolution with global coverage. This example
examines XCO2, the column-averaged mole fraction of atmospheric carbon
dioxide, from OCO-2. XCO2 is a scalar quantity reported in parts per
million (ppm) and is spatially and temporally referenced.

Most users of OCO-2 and OCO-3 XCO2 data work with the missions’ *lite*
products. These products are provided as global daily files at the
original satellite spatial resolution of approximately 2 km by 1 km. The
lite products include the XCO2 estimates, quality filtering information,
and other diagnostics. Data users who combine the global data with
geophysical models often choose to spatially aggregate these products,
and this example highlights some simple approaches to accomplish this
for quick visualization.

## Import packages

**One-time package configuration**

Earthdata authentication is handled using the `earthdatalogin` package.
Authentication for S3 credentials is handled using a development version
of the package. This installation should only need to be executed *once*
with the `pak` package.

``` r
library(pak)
pak::pkg_install("boettiger-lab/earthdatalogin")
```

Other packages are loaded for each R session.

``` r
suppressPackageStartupMessages({ 
    library(earthdatalogin)  
    library(httr)
    library(tidyr)
    library(dplyr)
    library(ggplot2)
    library(viridisLite)
    library(purrr)
    library(stringr)
    library(hdf5r)
    library(aws.s3)
})
```

## Authentication

Earthdatalogin authentication with username and password to obtain S3
credentials

``` r
edl_s3_token(daac = "https://data.gesdisc.earthdata.nasa.gov",prompt_for_netrc = FALSE)
```

## Search

The product archive can be searched efficiently using the Common
Metadata Repository (CMR) to identify products covering locations and
time periods of interest. This example searches a 3-day span of the
OCO-2 lite product archive.

- Make the CMR API request for OCO-2 LtCO2 products for Version 11.1r
- `httr` R package has similar capability to the Python `requests`
  library
- The `content` function parses the JSON response with the product
  information, producing a list of the granule information

``` r
cmr_url <- 'https://cmr.earthdata.nasa.gov/search/granules'

short_name <- 'OCO2_L2_Lite_FP'
start_time <- '2020-07-05T00:00:00Z'
end_time <- '2020-07-07T00:00:00Z'
time_string <- paste(start_time,end_time,sep=",")

cmr_query <- list(short_name=short_name,temporal=time_string,page_size="200",version="11.1r")
request_oco <- httr::GET(cmr_url,query=cmr_query)
# Confirm response in JSON format
httr::http_type(request_oco)
```

    ## [1] "application/json"

``` r
granules_v11 <- httr::content(request_oco,as="parsed") 
```

Now, the S3 URLs of the relevant products are extracted from the
`granules_v11` object. A custom function `href_s3` iterates through the
data links associated with each matched granule, and the desired S3 URLs
are retained. The result of this processing is a string array `s3_url`.

``` r
num_products <- length(granules_v11$feed$entry)
s3_url <- NULL

href_s3 <- function(link_index, ref_list) {
    # Extract the full link from a GES DISC JSON link list
    link_value <- ref_list[[link_index]]$href
    return(link_value)
}

for (i in seq(1,num_products)) {
    num_links <- length(granules_v11$feed$entry[[i]]$links)
    mapped_outputs <- unlist(map(seq(1,num_links), .f= href_s3, ref_list = granules_v11$feed$entry[[i]]$links))
    match_s3 <- stringr::str_detect(mapped_outputs,"s3://")
    s3_url <- c(s3_url,mapped_outputs[match_s3])
}

print(s3_url)
```

    ## [1] "s3://gesdisc-cumulus-prod-protected/OCO2_DATA/OCO2_L2_Lite_FP.11.1r/2020/oco2_LtCO2_200704_B11100Ar_230603215457s.nc4"
    ## [2] "s3://gesdisc-cumulus-prod-protected/OCO2_DATA/OCO2_L2_Lite_FP.11.1r/2020/oco2_LtCO2_200705_B11100Ar_230603215543s.nc4"
    ## [3] "s3://gesdisc-cumulus-prod-protected/OCO2_DATA/OCO2_L2_Lite_FP.11.1r/2020/oco2_LtCO2_200706_B11100Ar_230603215547s.nc4"
    ## [4] "s3://gesdisc-cumulus-prod-protected/OCO2_DATA/OCO2_L2_Lite_FP.11.1r/2020/oco2_LtCO2_200707_B11100Ar_230603215704s.nc4"

## S3 credential setup

- S3 direct access enabled by `aws.s3` package
- AWS region set with environmental variable

``` r
# Set AWS region
Sys.setenv("AWS_DEFAULT_REGION" = "us-west-2")
```

## Open granule

- A single OCO-2 lite product is accessed using the tools from `aws.s3`
  and `hdf5r`. The `ls` function lists the top-level variables (and
  groups) in the product.
- In addition to geolocation information and the XCO2 data, the
  `sounding_id` is extracted. OCO-2/3 identify each observation with a
  unique sounding ID, which also encodes a timestamp of the observation.
- Observation quality information is available with the
  `xco2_quality_flag`, where values of 0 are assigned to the best
  quality observations. See also [OCO-2/3 data user’s
  guide](https://docserver.gesdisc.eosdis.nasa.gov/public/project/OCO/OCO2_V11.2_OCO3_V11_L2_Data_Users_Guide_250304.pdf)

``` r
# S3 access dataset
file_object <- aws.s3::s3read_using(FUN = hdf5r::H5File$new, object = s3_url[1])
file_object$ls(recursive=FALSE)$name
```

    ##  [1] "xco2_apriori"           "vertex_latitude"        "file_index"            
    ##  [4] "vertices"               "xco2_qf_simple_bitflag" "vertex_longitude"      
    ##  [7] "pressure_levels"        "xco2"                   "source_files"          
    ## [10] "time"                   "pressure_weight"        "Preprocessors"         
    ## [13] "solar_zenith_angle"     "longitude"              "xco2_x2019"            
    ## [16] "xco2_qf_bitflag"        "latitude"               "sensor_zenith_angle"   
    ## [19] "levels"                 "Meteorology"            "xco2_quality_flag"     
    ## [22] "sounding_id"            "xco2_averaging_kernel"  "Auxiliary"             
    ## [25] "bands"                  "date"                   "Retrieval"             
    ## [28] "Sounding"               "xco2_uncertainty"       "co2_profile_apriori"   
    ## [31] "epoch_dimension"

``` r
lat <- file_object[['latitude']][]
lon <- file_object[['longitude']][]
quality_flag_xco2 <- file_object[['xco2_quality_flag']][]
xco2 <- file_object[['xco2']][]
sounding_id <- file_object[['sounding_id']][]
file_object$close_all()

table(quality_flag_xco2)
```

    ## quality_flag_xco2
    ##     0     1 
    ## 91004 67297

- These OCO-2 variables are assembled into a data frame known as a
  *tibble*, which is a modernized data frame object that is part of the
  [tidyverse](https://tibble.tidyverse.org/)
- A new variable `Sounding_10sec` is created, which truncates the
  sounding ID into a 10-second interval for later grouping and
  aggregation. Due to the OCO-2 sampling frequency, up to 240
  observations could be grouped into the same 10-second interval.
- The data are subset to include only high-quality observations.

``` r
xco2_tibble <- tibble(SoundingID=as.vector(sounding_id), Latitude=as.vector(lat), Longitude=as.vector(lon),
                XCO2=as.vector(xco2), V11QFlag=as.vector(quality_flag_xco2))
xco2_tibble <- xco2_tibble %>% dplyr::mutate(Sounding_10sec = floor(SoundingID / 1.0e3))
xco2_tibble <- xco2_tibble %>% dplyr::filter(V11QFlag == 0)

print(head(xco2_tibble))
```

    ## # A tibble: 6 × 6
    ##   SoundingID Latitude Longitude  XCO2 V11QFlag Sounding_10sec
    ##        <dbl>    <dbl>     <dbl> <dbl>    <int>          <dbl>
    ## 1    2.02e15    -38.3     -152.  412.        0  2020070400004
    ## 2    2.02e15    -36.7     -153.  413.        0  2020070400011
    ## 3    2.02e15    -34.6     -153.  411.        0  2020070400015
    ## 4    2.02e15    -34.6     -153.  410.        0  2020070400015
    ## 5    2.02e15    -34.6     -153.  412.        0  2020070400015
    ## 6    2.02e15    -34.6     -153.  411.        0  2020070400015

## Group and Summarize

Construct a summary tibble that provides the median longitude, latitude,
and XCO2 of observations within each 10-second block. Since the same
summary function will be applied to multiple variables, this task can
make use of the data wrangling functions in the `tidyr` package.
Specifically, the [longer/wider reshaping
procedure](https://jhudatascience.org/tidyversecourse/wrangle-data.html)
is applied

- Reshape tibble to *longer* format with variable and value columns
  using `pivot_longer`
- Group the result by variable and 10-second interval and compute
  summaries and retain only groups with over 30 observations
- Reshape grouped data back to *wider* format for additional analysis and
  plotting using `pivot_wider`

``` r
# Use the pivot_longer, pivot_wider approach
xco2_tibble_longer <- xco2_tibble %>% tidyr::pivot_longer(cols = c("XCO2","Latitude","Longitude"), names_to = "GeoVar",
                                      names_prefix="", values_to = "value")
xco2_group <- xco2_tibble_longer %>% dplyr::group_by(Sounding_10sec,GeoVar) %>% 
                                     dplyr::summarise(Med = median(value,na.rm=TRUE), NumSamples=n()) %>% dplyr::ungroup()
```

    ## `summarise()` has grouped output by 'Sounding_10sec'. You can override using
    ## the `.groups` argument.

``` r
xco2_group <- xco2_group %>% dplyr::filter(NumSamples > 30)
xco2_group_wider <- xco2_group %>% tidyr::pivot_wider(id_cols = c("Sounding_10sec"), names_from=c("GeoVar"), 
                                               names_prefix=c("Med_"), values_from="Med")
```

## Plot

- Aggregated XCO2 data are mapped with the `ggplot2` package and a
  custom plot theme `theme_mat`
- Tibble is converted to a `sf` object to encode geospatial coordinate
  information for mapping

``` r
# Map results
suppressPackageStartupMessages({ 
    library("rnaturalearth")
    library("rnaturalearthdata")
    library(sf)
})


theme_mat = ggplot2::theme_bw() 
theme_mat$axis.title.x$size = 11
theme_mat$axis.text.x$size = 10
theme_mat$axis.title.y$size = 11
theme_mat$axis.text.y$size = 10
theme_mat$plot.title$size = 12
theme_mat$legend.position = "right"
theme_mat$plot.title$hjust = 0.5
theme_mat$strip.text$size = 11
theme_mat$legend.text$size = 10
theme_mat$legend.title$size = 11

xco2_merged <- as.data.frame(xco2_group_wider)
xco2_merged_sf <- sf::st_as_sf(xco2_merged, coords = c("Med_Longitude", "Med_Latitude"), 
                               crs = 4326, agr = "constant")

world <- rnaturalearth::ne_coastline(scale = "medium", returnclass = "sf")
p9 = viridisLite::viridis(9)

xco2_map <- ggplot(data = world) + geom_sf(color = '#777777') + 
  geom_sf(data = xco2_merged_sf,aes(color=Med_XCO2), size=0.8) +
  scale_color_gradientn("XCO2 [ppm]",colors=p9)  + 
  theme_mat + ggtitle("OCO-2 10 sec Median of XCO2")  

xco2_map
```

![](Rfigures/map10s-1.png)<!-- -->