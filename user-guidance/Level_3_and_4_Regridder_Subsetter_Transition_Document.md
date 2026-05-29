# GES DISC Level 3 & 4 Regridder and Subsetter Replaced with Enterprise Level 3 & 4 Subsetter

In line with NASA’s goal of unification with Enterprise services, the GES DISC Level 3 & 4 Regridder and Subsetter Service is being transitioned into an integrated Enterprise-powered Level 3 and 4 Subsetting Service. This transition will take place no earlier than July 31st, 2026. This document outlines the capabilities of the new Level 3 and 4 Subsetting Service, including links to relevant How-Tos.

## Capabilities
**The new Level 3 & 4 Subsetting Service provides the following capabilities:**

- Subset a bounding box
- Subset a point-radius
- Subset by a shape or polygon
- Subset by variables of interest
- Subset by dimension (supported by Harmony API only)
- Temporal subsetting
- Preserve original file structure, metadata, dimensions, data types, and attributes
  
**The following capabilities are currently not supported in the new subsetting service:**

- Subset by recurring time of day
- Regridding
- Output as GeoTIFF/COG
- Statistics (mean/max/min)
- Select custom dimension range

## How Can Users Access the New Level 3 & 4 Subsetting Service?
The new Level 3 & 4 Subsetting Service uses the Harmony Application Programming Interface (API) in the [Earthdata Cloud](https://www.earthdata.nasa.gov/about/earthdata-cloud-evolution) to search and perform subsetting. Users will need to be logged into their Earthdata Login account to access the new Level 3 & 4 Subsetting Service. The following three methods can be used to perform subsetting:

### 1. *Level 3 & 4 Subsetting Service User Interface:*
For select Level 3 & 4 collections, a user interface is available on their associated dataset landing page beneath the "Data Access" panel. This method requires an [Earthdata Login](https://urs.earthdata.nasa.gov/documentation/what_do_i_need_to_know) account. To learn more about how to use the subsetter user interface, please visit this how-to: [How to Subset Level 2 Data with the Earthdata Enterprise Subsetter](http://disc.gsfc.nasa.gov/information/howto?title=How%20to%20Subset%20Level%202%20Data%20with%20the%20Earthdata%20Enterprise%20Subsetter), which outlines how to use the subsetter for Level 2 data, but these steps can be applied to  Level 3 and 4 collections. 

### 2. *Harmony API*

The Harmony API handles various subsetting operations on select Level 3 & 4 data, running from the [Earthdata Cloud](https://www.earthdata.nasa.gov/about/earthdata-cloud-evolution) and accessible to anyone with an [Earthdata Login](https://urs.earthdata.nasa.gov/documentation/what_do_i_need_to_know) account. It can be accessed by querying URLs that utilize Open Geospatial Consortium (OGC) syntax, or through the Python [Harmony-py](http://harmony-py.readthedocs.io/en/main/) package, as demonstrated in this Jupyter Notebook tutorial: [How to Access Level 2, 3, and 4 Data Using Python](https://disc.gsfc.nasa.gov/information/howto?title=How%20to%20Subset%20and%20Download%20Level%202,%203,%20and%204%20GES%20DISC%20Data%20Using%20Python).

### 3. *Earthdata Search*

The Earthdata Search website, managed by NASA ESDIS provides a user interface for searching and subsetting collections from all NASA Earthdata Data Access and Archive Centers (DAACs), including the ability to subset select Level 3 & 4 GES DISC collections using the Harmony API Interface. This method requires an Earthdata Login account. To learn more about how to use the subsetter user interface, please visit this how-to:  [How do I find data using Earthdata Search?](https://nasa-openscapes.github.io/earthdata-cloud-cookbook/how-tos/find-data/earthdata_search.html)

## Guides and Resources

- For guides on using the Dataset Landing Page subsetter user interface:
    - [How to Subset Level 2 Data with the Earthdata Enterprise Subsetter ](http://disc.gsfc.nasa.gov/information/howto?title=How%20to%20Subset%20Level%202%20Data%20with%20the%20Earthdata%20Enterprise%20Subsetter)
- For guides on using the Earthdata Search website:
    - [How do I find data using Earthdata Search?](https://nasa-openscapes.github.io/earthdata-cloud-cookbook/how-tos/find-data/earthdata_search.html)
    - [Video: Using Harmony Tools in Earthdata Search](https://www.youtube.com/watch?v=nrs2wTbcp-M)
- For guides on accessing the service using programming tools:
    - [Python: How to Subset and Download Level 2, 3 and 4 Data using Harmony-py](http://disc.gsfc.nasa.gov/information/howto?title=How%20to%20Subset%20and%20Download%20Level%202,%203,%20and%204%20GES%20DISC%20Data%20Using%20Python)
      
## Questions?
Please reach out on the [Earthdata Forum](https://forum.earthdata.nasa.gov/viewforum.php?f=7) or contact the GES DISC Help Desk (via email at gsfc-dl-help-disc@mail.nasa.gov) for any questions about cloud data and services. Other contact methods are listed on our Contact Us page.
