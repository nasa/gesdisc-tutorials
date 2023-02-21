# GES DISC Cloud Tutorials

A place to find cloud-relevant tutorials on how to use GES DISC tools, services, and data.

Most tutorials in this repository take the form of Python notebooks. Jupyter is a very popular version of Python notebooks, and is used extensively by the GES DISC team.

## What You Need To Know About the Earthdata Cloud

Keep up with current cloud options by visiting the [Cloud Migration page](https://disc.gsfc.nasa.gov/information/documents?title=Migrating%20to%20the%20Cloud) on the GES DISC website. Here you will find [basic information on the transition to the cloud](https://disc.gsfc.nasa.gov/information/documents?title=Migrating%20to%20the%20Cloud#introduction), [frequently asked questions](https://disc.gsfc.nasa.gov/information/documents?title=Migrating%20to%20the%20Cloud#faq), [How-Tos](https://disc.gsfc.nasa.gov/information/documents?title=Migrating%20to%20the%20Cloud#how-to) and [definitions](https://disc.gsfc.nasa.gov/information/glossary?keywords=%22Earthdata%20Cloud%22&page=1). 

Some tutorials can be run locally to take advantage of cloud archives or cloud-based tools, while other tutorials must be executed using cloud resources provided by AWS in a particular region, us-west-2. These tutorials contain a salmon-colored banner indicating this requirement.

![](../images/us-west-2-banner.png "Banner used in notebooks that require compute resources in the AWS us-west-2 region.")



| Notebook  | Summary | Services and Tools |
| ------------- |-------------|:-------------:|
|[How to Create and Store Earthdata Login Credentials Using Pthon](notebooks/How_to_Create_and_Store_Earthdata_Login_Credentials_Using_Python.ipynb) | This notebook demonstrates how to generate and store your Earthdata Login credentials in a .netrc file. | |
|[How to Directly Access MERRA-2 Data from an S3 Bucket with Python](notebooks/How_to_Directly_Access_MERRA-2_Data_from_an_S3_Bucket.ipynb) | This notebook demonstrates how to access and plot a Modern-Era Retrospective analysis for Research and Applications (MERRA-2) M2T1NXSLV.5.12.4 file hosted via an Amazon S3 bucket. It demonstrates how to access an S3 bucket with the S3FS library and then plot sea-level pressure contours of a single file with Cartopy and Matplotlib.| |
|[How to Obtain a List of S3 URLs for a GES DISC Collection using the CMR API](/notebooks/How_to_Directly_Access_MERRA-2_Data_from_an_S3_Bucket.ipynb)| This notebook demonstrates how to obtain a list of S3 URLs for desired cloud-hosted GES DISC granules using the Commmon Metadata Repository (CMR) API. |  |
|[How to Retrieve Temporary S3 Credentials for the GES DISC Cloud Archive](notebooks/How_to_Retrieve_Temporary_S3_Credentials_for_the_GES_DISC_Cloud_Archive.ipynb) | This notebook demonstrates how to retrieve GES DISC S3 credentials by using a previously generated netrc file.  | |
