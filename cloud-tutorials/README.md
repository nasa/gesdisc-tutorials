# GES DISC Cloud Tutorials

A place to find cloud-relevant tutorials on how to use GES DISC tools, services, and data.

Most tutorials in this repository take the form of Python notebooks. Jupyter is a very popular version of Python notebooks, and is used extensively by the GES DISC team.

## What You Need To Know About the Earthdata Cloud

Keep up with current cloud options by visiting the [Cloud Migration page](https://disc.gsfc.nasa.gov/information/documents?title=Migrating%20to%20the%20Cloud) on the GES DISC website. Here you will find [basic information on the transition to the cloud](https://disc.gsfc.nasa.gov/information/documents?title=Migrating%20to%20the%20Cloud#introduction), [frequently asked questions](https://disc.gsfc.nasa.gov/information/documents?title=Migrating%20to%20the%20Cloud#faq), [How-Tos](https://disc.gsfc.nasa.gov/information/documents?title=Migrating%20to%20the%20Cloud#how-to) and [definitions](https://disc.gsfc.nasa.gov/information/glossary?keywords=%22Earthdata%20Cloud%22&page=1). 

Some tutorials can be run locally to take advantage of cloud archives or cloud-based tools, while other tutorials must be executed using cloud resources provided by AWS in a particular region, us-west-2. These tutorials contain a salmon-colored banner indicating this requirement.

![](../images/us-west-2-banner.png "Banner used in notebooks that require compute resources in the AWS us-west-2 region.")



| Notebook  | Summary | Services, Tools, Data Types | Actions |
| ------------- |-------------|:-------------:|:-------------:|
|[How to Create and Store Earthdata Login Credentials Using Python](notebooks/How_to_Create_and_Store_Earthdata_Login_Credentials_Using_Python.ipynb) | This notebook demonstrates how to generate and store your Earthdata Login credentials in a .netrc file. | Python | Access |
|[How to Directly Access MERRA-2 Data from an S3 Bucket with Python](notebooks/How_to_Directly_Access_MERRA-2_Data_from_an_S3_Bucket.ipynb) | This notebook demonstrates how to access and plot a Modern-Era Retrospective analysis for Research and Applications (MERRA-2) M2T1NXSLV.5.12.4 file hosted via an Amazon S3 bucket. It demonstrates how to access an S3 bucket with the S3FS library and then plot sea-level pressure contours of a single file with Cartopy and Matplotlib.| Python, Direct S3 Access | Access, Subset, Plot |
|[How to Directly Access GPM IMERG Data from an S3 Bucket with Python](notebooks/How_to_Directly_Access_GPM_IMERG_Data_from_an_S3_Bucket.ipynb) | This notebook demonstrates how to search, access, and plot an Integrated Multi-satellitE Retrievals for GPM file hosted via an Amazon S3 bucket. It demonstrates how to search for a granule by DOI and access it directly from an S3 bucket using the earthaccess and Xarray libraries, before plotting precipitation for that granule using Cartopy and Matplotlib.| Python, Direct S3 Access | Search, Access, Subset, Plot |
|[How to Obtain a List of S3 URLs for a GES DISC Collection using the CMR API](notebooks/How_to_Obtain_a_List_of_S3_URLs_for_a_GES_DISC_Collection_Using_the_CMR_API.ipynb)| This notebook demonstrates how to obtain a list of S3 URLs for desired cloud-hosted GES DISC granules using the Commmon Metadata Repository (CMR) API. | Python, CMR | Search |
|[How to Retrieve Temporary S3 Credentials for the GES DISC Cloud Archive](notebooks/How_to_Retrieve_Temporary_S3_Credentials_for_the_GES_DISC_Cloud_Archive.ipynb) | This notebook demonstrates how to retrieve GES DISC S3 credentials by using a previously generated netrc file.  | Python | Access |
|[How to Perform Cross-DAAC S3 Bucket Access Using Python](notebooks/How_to_Perform_Cross-DAAC_S3_Bucket_Access_Using_Python.ipynb) | This notebook demonstrates how to access cloud-hosted Earthdata granules from S3 buckets using the CMR API and Python, from two different DAACs (GES DISC and PO.DAAC).  | Python, CMR, Direct S3 Access | Search, Access, Subset, Compute, Plot |
|[How to Search and Access Giovanni Variable Zarr Stores](notebooks/How_to_Search_and_Load_Zarr_Stores.ipynb) | This notebook demonstrates how to search for and access Giovanni Cache Zarr Stores inside the AWS us-west-2 region.  | Python, CMR, Direct S3 Access | Search, Access, Subset, Compute, Plot |
|[How to Access and Analyze GPM Data from Giovanni Variable Zarr Stores](notebooks/How_to_Access_and_Analyze_GPM_Zarr_Store.ipynb) | This notebook demonstrates how to access GPM data from the Giovanni Cache Zarr Store inside the AWS us-west-2 region for performing analysis.  | Python, CMR, Direct S3 Access | Search, Access, Subset, Compute, Plot |
