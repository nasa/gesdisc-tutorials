#!/usr/bin/env python3
"""
- This script performs collocation analysis between Level 2 (L2) TROPOMI satellite aerosol optical depth (AOD) measurements and ground-based AERONET observations for validation and comparison purposes. 

- The tool matches satellite pixels within 25 km of each AERONET site with AERONET observations averaged within ±0.5 hours of the satellite overpass time (these matching criteria are adjustable.) 

- Users can add additional satellite datasets as needed by following the TROPOMI implementation example. 

### Author: Xiaohua Pan
- Email: xiaohua.pan@nasa.gov
- Organization: NASA Goddard Space Flight Center
- Date Created: 09/02/2025
- Last Modified: 11/17/2025

### Features:
- Downloads TROPOMI data from NASA Earthdata
- Retrieves AERONET data from web services
- Performs spatial and temporal matching with user-defined criteria
- Generates visualization and statistical analysis outputs

### Data Collocation Process 
#### Step 1: Establish Collocation Criteria (User-Adjustable) 
- Default Spatial Criteria: Satellite pixels within a 25 km radius of each AERONET site
- Default Temporal Criteria: AERONET observations averaged within ±0.5 hours of satellite overpass time
- Sample Use Case: AOD analysis during Los Angeles wildfire events (January 8-12, 2025)

  #### Implementation Note:
- Initial analyses should utilize small spatial domains (approximately 2° × 2°) to minimize processing time
- Short temporal periods (5 days or less) are recommended for initial testing

#### Step 2: Data Acquisition
- AERONET Data: Retrieve AOD observations from sites within the defined spatial domain and temporal range
- TROPOMI Data: Access L2 AOD swath data covering the corresponding spatial and temporal parameters

#### Step 3: Spatial Collocation
- Identify TROPOMI L2 AOD pixels falling within 25 km of each AERONET site
- Calculate statistical metrics (mean, standard deviation) of AOD values from collocated pixels
- Define satellite overpass time as the mean acquisition time of collocated pixels within 25 km of each AERONET site

#### Step 4: Temporal Collocation
- Match satellite observations with AERONET data averaged within ±0.5 hours of satellite overpass time
- Create paired satellite-ground observation datasets using datetime as satellite passing time

#### Step 5: Visualization and Analysis
- Collocation Map: Generate spatial representation of successfully collocated AERONET sites
- Time Series Analysis: Produce temporal plots of collocated AERONET and TROPOMI L2 AOD measurements at each site
- Correlation Analysis: Create scatter plots comparing collocated measurements for individual sites and aggregate data across all locationsStep 

### Input data:
- AERONET: version 3, level 1.5 or level 2 (https://aeronet.gsfc.nasa.gov/print_web_data_help_v3_new.html)
-- Temporal coverage: vary from site to site
-- Temporal resolution: every 3-15 minutes
-- Spatial coverage: global

- TROPOMAER: version 1, level 2 (https://doi.org/10.5067/MEASURES/AER/DATA204)
-- Temporal coverage: 2018-04-30 to present
-- Temporal resolution: once per day
-- Spatial coverage: global
-- Spatial resolution: 7.5km x 3km 

### Output:
All outputs are saved in the output_dir
- /plots: Map of AERONET sites with valid data; scatter plots (aggregate and per-site) and time series (per-site)
- /data: Collocation results exported in CSV format
- /logs: initialization.log 

### Requirements:
- Earthdata credentials (https://disc.gsfc.nasa.gov/earthdata-login) 
- Earthdata prerequisite Files (https://disc.gsfc.nasa.gov/information/howto?title=How%20to%20Generate%20Earthdata%20Prerequisite%20Files)
- Python packages: earthaccess, pandas, numpy, matplotlib, requests, BeautifulSoup4, netCDF4

### Tested Configurations:
This notebook was tested succussfully in the configurations below
- Mac OS, Python 3.11.0, NetCDF4 1.6.2 
- Red Hat Enterprise Linux 8.10, Python 3.12.2, NetCDF4 1.6.2
- Google colab 
### Caveats:
- It is possible that no AERONET sites or TROPOMI swaths fall within the specified search domain.
- Collocated observations may not be found for the selected criteria.
- AERONET data downloads may fail due to server limitations (timeout errors—retry later)
- Downloaded TROPOMI files may be corrupted; verify before processing
- Use a text viewer (e.g., vim, not Microsoft Excel) to view downloaded CSV AERONET data in ./data to ensure dates are in DD:MM:YYYY format

### Version History:
- v1  - [09/12/2025] - [Xiaohua Pan]
    - Initial stakeholder demonstration
    
- v2  - [09/25/2025] - [Xiaohua Pan]
    - Match the satellite with the nearest AERONET observation within ±0.5 hours of satellite passing time.
    
- v3 - [10/03/2025] - [Xiaohua Pan]
    - Bin AERONET measurements into 1-hour intervals and collocates the satellite with the closest bin  
    
- v4 - [10/30/2025] - [Xiaohua Pan]
    - Match the satellite with the AERONET observation averaged within ±0.5 hours of satellite overpass time

- v5 - [11/17/2025] - [Xiaohua Pan]
    - Match the satellite with the AERONET observation averaged within ±0.5 hours of satellite overpass time
    - Fixed some collocation issues found in v4. 

### Acknowledgments:
- ChatGSFC AI assistant used for code improvement and debugging
- NASA GES DISC colleagues, especially Thomos Hearty for validating and Christopher Battisto for publishing this notebook. 
- TROPOMI aerosol data provided by ESA/NASA MEaSUREs program
- AERONET data courtesy of NASA/GSFC AERONET team
    
### Disclaimer:
This software is provided "as is" without any warranty of any kind, either 
expressed, implied, or statutory, including, but not limited to, any warranty 
that the software will conform to specifications, any implied warranties of 
merchantability, fitness for a particular purpose, and freedom from infringement.

"""

import os
import io
#import sys
import glob
import shutil
import numpy as np
import pandas as pd
import netCDF4 as nc
#from bs4 import BeautifulSoup as bso
from datetime import datetime, timedelta
import requests
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
#from scipy.spatial.distance import cdist
import warnings
import earthaccess
warnings.filterwarnings('ignore')

class TropomiAeronetCollocation:
    def __init__(self, spatial_radius_km=25.0, temporal_window_hours=0.5, output_dir='./collocation_results', 
                 target_wavelength=500, start_date=None, end_date=None, bounding_box=None, 
                 aeronet_data_type='AOD15', aeronet_data_format=10, tropomi_data_shortname='TROPOMAER', tropomi_data_version=1):
        """
        Initialize the TROPOMI-AERONET collocation class
        
        Parameters:
        -----------
        spatial_radius_km : float
            Maximum distance in kilometers for spatial collocation
        temporal_window_hours : float
            Maximum time difference in hours for temporal collocation
        output_dir : str
            Directory to save results and plots
        target_wavelength : int
            Target wavelength for AERONET AOD comparison (e.g., 380, 500, 675 nm)
        start_date : datetime
            Start date for analysis period
        end_date : datetime
            End date for analysis period
        bounding_box : dict, required
            Dictionary with keys 'min_lon', 'max_lon', 'min_lat', 'max_lat' 
            defining the region to search for AERONET sites.
            Longitude values must be in range [-180, 180].
            Latitude values must be in range [-90, 90].
            Example: {'min_lon': -120, 'max_lon': -115, 'min_lat': 32, 'max_lat': 35}
        aeronet_data_type : str
            AERONET data type (e.g., 'AOD20' = Level 2.0; 'AOD15' = Level 1.5)
        aeronet_data_format : int
            AERONET data format (e.g., 0 = all points; 10 = daily average)
        tropomi_data_shortname: str
            tropomi AOD product's shortname, such as 'TROPOMAER'
        tropomi_data_shortname: str
                tropomi AOD product's version number, such as '1'
            
        """
        self.spatial_radius_km = spatial_radius_km
        self.temporal_window_hours = temporal_window_hours
        self.output_dir = output_dir
        self.target_wavelength = target_wavelength
        self.start_date = start_date
        self.end_date = end_date
        self.bounding_box = bounding_box
        self.aeronet_data_type = aeronet_data_type
        self.aeronet_data_format = aeronet_data_format
        self.tropomi_data_shortname = tropomi_data_shortname
        self.tropomi_data_version = tropomi_data_version
        
        # Remove output directory completely if it exists
        if os.path.exists(self.output_dir): shutil.rmtree(self.output_dir)

        # Create fresh output directory and subdirectories
        os.makedirs(self.output_dir, exist_ok=True)
        # Create subdirectories for organizing output
        os.makedirs(os.path.join(self.output_dir, "plots"), exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "data"), exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "logs"), exist_ok=True)
        
        # Initialize empty dictionary for AERONET sites
        self.aeronet_sites = {}
        
      # Print initialization summary
        print("="*60)
        print("TROPOMI-AERONET Collocation Analysis Initialized")
        print("="*60)
        print(f"  Spatial radius: {self.spatial_radius_km} km")
        print(f"  Temporal window: ±{self.temporal_window_hours} hours")
        print(f"  Target wavelength: {self.target_wavelength} nm")
        print(f"  Output directory: {self.output_dir}")
        
        if self.start_date and self.end_date:
            print(f"  Analysis period: {self.start_date.strftime('%Y-%m-%d')} to {self.end_date.strftime('%Y-%m-%d')}")
        
        if self.bounding_box:
            print(f"  Region of interest:")
            print(f"    Longitude: {self.bounding_box['min_lon']} to {self.bounding_box['max_lon']}")
            print(f"    Latitude: {self.bounding_box['min_lat']} to {self.bounding_box['max_lat']}")
        
        print(f"  AERONET data type: {self.aeronet_data_type}")
        print(f"  AERONET data format: {self.aeronet_data_format}")
        
        
        print(f"  TROPOMI product's shortname': {self.tropomi_data_shortname}")
        print(f"  TROPOMI product's version': {self.tropomi_data_version}")
        
        
        # Search for AERONET sites in the specified region
        self.aeronet_sites = self.search_aeronet_sites()
        
        # Warning if no sites found
        if not self.aeronet_sites:
            print("="*60)
            print("WARNING: No AERONET sites found in the specified bounding box.")
            print("Consider the following options:")
            print("  1. Expand your bounding box")
            print("  2. Check a different region")
            print("  3. Verify AERONET site availability at:")
            print("     https://aeronet.gsfc.nasa.gov/cgi-bin/draw_map_display_aod_v3")
            print("="*60)
            return
        
        # For backward compatibility with existing code
        #self.la_aeronet_sites = self.aeronet_sites
        
        # Print configuration summary
        site_count = len(self.aeronet_sites)
        print(f"\nSuccessfully initialized TROPOMI-AERONET collocation system")
        print(f"  Found {site_count} AERONET sites in region")
        print(f"  Spatial radius: {self.spatial_radius_km} km")
        print(f"  Temporal window: ±{self.temporal_window_hours} hours")
        #print(f"  Analysis period: {self.analysis_days} days")
        print(f"  Output directory: {self.output_dir}")
        
        # Print the sites being used
        if site_count > 0:
            print("\nAERONET sites in selected region:")
            for site_name, coords in self.aeronet_sites.items():
                elev_str = f", elev={coords['elevation']} m" if 'elevation' in coords and coords['elevation'] is not None else ""
                print(f"  - {site_name}: lat={coords['lat']:.4f}°, lon={coords['lon']:.4f}°{elev_str}")
        

       # Search and download TROPOMI data in the specified region
        self.tropomi_data_dir = self.download_tropomi_data()
        
        # Add this fallback handling
        if self.tropomi_data_dir is None:
            # Fallback to a default directory if download failed
            self.tropomi_data_dir = os.path.join(self.output_dir, 'data/TROPOMI_fallback')
            print(f"Download failed, using fallback directory: {self.tropomi_data_dir}")
       
        # TROPOMI data files (modify path as needed)
        #tropomi_data_dir = "/Users/xpan2/Library/CloudStorage/OneDrive-NASA/work2025/Python/JNbook/MAPSS/data/TROPOMAER.1"
        
        if not os.path.exists(self.tropomi_data_dir):
            print(f"\nCreating example directory structure at {self.tropomi_data_dir}")
            print("Please place your TROPOMI L2 aerosol netCDF files in this directory")
            os.makedirs(self.tropomi_data_dir, exist_ok=True)
        
        # Find all TROPOMI files in the directory
        self.tropomi_files = sorted(glob.glob(f"{self.tropomi_data_dir}/*.nc*"))
        
        if len(self.tropomi_files) == 0:
            print(f"\nNo TROPOMI netCDF files found in {self.tropomi_data_dir}")
            print("Please download TROPOMI L2 aerosol files and place them in this directory")
            return
        
        print(f"\nFound {len(self.tropomi_files)} TROPOMI files in {self.tropomi_data_dir}")
            
        # Log initialization
        self.log_initialization()
        

    
    def log_initialization(self):
        """Log initialization parameters to a file for record-keeping"""
        log_file = os.path.join(self.output_dir, "logs", "initialization.log")
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        with open(log_file, 'w') as f:
            f.write(f"TROPOMI-AERONET Collocation Initialization\n")
            f.write(f"Timestamp: {timestamp}\n")
            f.write(f"--------------------------------------\n")
            f.write(f"Output directory: {self.output_dir}\n")
            
            f.write(f"\nCollocation Parameters:\n")
            f.write(f"  - Spatial radius of an AERONET site: {self.spatial_radius_km} km\n")
            f.write(f"  - Temporal window: ±{self.temporal_window_hours} hours - then match the closest AERONET data to the satelite passing time \n")
            f.write(f"  - Target wavelength: {self.target_wavelength} nm\n")
            
            f.write(f"\nAnalysis Period:\n")
            if self.start_date and self.end_date:
                f.write(f"  - Start date: {self.start_date.strftime('%Y-%m-%d')}\n")
                f.write(f"  - End date: {self.end_date.strftime('%Y-%m-%d')}\n")
                days_diff = (self.end_date - self.start_date).days
                f.write(f"  - Duration: {days_diff} days\n")
            else:
                f.write(f"  - Start date: Not specified\n")
                f.write(f"  - End date: Not specified\n")
            
            f.write(f"\nRegion of Interest:\n")
            if self.bounding_box:
                f.write(f"  - Longitude range: {self.bounding_box['min_lon']} to {self.bounding_box['max_lon']}\n")
                f.write(f"  - Latitude range: {self.bounding_box['min_lat']} to {self.bounding_box['max_lat']}\n")
            else:
                f.write(f"  - Bounding box: Not specified\n")
            
            f.write(f"\nAERONET Configuration:\n")
            f.write(f"  - Data type: {self.aeronet_data_type}\n")
            f.write(f"  - Data format: {self.aeronet_data_format}\n")
            
            f.write(f"\nFound AERONET sites ({len(self.aeronet_sites)}):\n")
            f.write("Note: not all sites have requested data. Check the Collocation Summary print to the screen and the plots")
            if self.aeronet_sites:
                for site_name, coords in self.aeronet_sites.items():
                    elev_str = f", elev={coords['elevation']} m" if 'elevation' in coords and coords['elevation'] is not None else ""
                    f.write(f"  - {site_name}: lat={coords['lat']:.4f}°, lon={coords['lon']:.4f}°{elev_str}\n")
            else:
                f.write(f"  - No sites found in specified region\n")
                
            f.write(f"\nTROPOMI Configuration:\n")
            f.write(f"  - Data shortname: {self.tropomi_data_shortname}\n")
            f.write(f"  - Data version: {self.tropomi_data_version}\n")
            
            if self.tropomi_files:
                f.write(f"\nFound {len(self.tropomi_files)} TROPOMI files \n")
            else:
                f.write(f"  - No TROPOMI files found in specified region\n")
                  
        
        print(f"Initialization details logged to: {log_file}")
    
    def haversine_distance(self, lat1, lon1, lat2, lon2):
        """
        Calculate the great circle distance between two points in km
        """
        R = 6371.0  # Earth's radius in km
        
    # Convert to radians
        lat1, lon1 = np.radians(lat1), np.radians(lon1)
        #lat2_array, lon2_array = np.radians(lat2_array), np.radians(lon2_array)
        lat2_array, lon2_array = np.radians(lat2), np.radians(lon2)
       
       # Vectorized calculation
        dlat = lat2_array - lat1
        dlon = lon2_array - lon1
       
        a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2_array) * np.sin(dlon/2)**2
        c = 2 * np.arcsin(np.sqrt(a))
        
        return R * c
     
    def search_aeronet_sites(self):
        """
        Search for AERONET sites within a specified bounding box
        
        Parameters:
        -----------
        min_lon, max_lon : float
            Minimum and maximum longitude (degrees), in range [-180, 180]
        min_lat, max_lat : float
            Minimum and maximum latitude (degrees), in range [-90, 90]
            
        Returns:
        --------
        dict
            Dictionary of sites within the bounding box, with site names as keys
            and dictionaries with 'lat', 'lon', 'elevation' as values
        """
        print(f"Downloading AERONET site list...")
        
        # URL for AERONET site locations file
        #aeronet_sites_url = "https://aeronet.gsfc.nasa.gov/aeronet_locations_v3.txt"   # all sites 
        aeronet_sites_url = "https://aeronet.gsfc.nasa.gov/aeronet_locations_v3921.txt"  # only active sites
        
        min_lon = self.bounding_box['min_lon'] 
        max_lon = self.bounding_box['max_lon']
        min_lat = self.bounding_box['min_lat']
        max_lat = self.bounding_box['max_lat']
        
        try:
            # Download the site list
            response = requests.get(aeronet_sites_url, timeout=30)
            
            if response.status_code != 200:
                print(f"Failed to download AERONET site list. Status code: {response.status_code}")
                return {}
                
            # Parse the site information
            sites_within_bbox = {}
            lines = response.text.split('\n')
            
            # Skip header lines and metadata lines
            data_lines = []
            header_found = False
            
            for line in lines:
                line = line.strip()
                if not line:
                    continue
                    
                # Skip metadata lines that contain '=' or start with '#'
                if '=' in line or line.startswith('#'):
                    continue
                    
                # Skip the header line with column names
                if 'Site_Name' in line and 'Longitude' in line and 'Latitude' in line:
                    header_found = True
                    continue
                    
                # Only process data lines after we've found the header
                if header_found:
                    data_lines.append(line)
            
            # Handle special case where min_lon > max_lon (crossing 180/-180 boundary)
            crosses_date_line = min_lon > max_lon
            
            for line in data_lines:
                try:
                    # AERONET locations format is now:
                    # Site_Name,Longitude(decimal_degrees),Latitude(decimal_degrees),Elevation(meters)
                    parts = [p.strip() for p in line.split(',')]
                    
                    if len(parts) < 3:
                        continue
                    
                    site_name = parts[0].strip()
                    lon = float(parts[1])  # Note: longitude comes before latitude in the new format
                    lat = float(parts[2])
                    
                    # Get elevation if available
                    elevation = None
                    if len(parts) >= 4 and parts[3]:
                        try:
                            elevation = float(parts[3])
                        except ValueError:
                            elevation = None
                    
                    # Check if site is within the bounding box
                    lat_in_range = min_lat <= lat <= max_lat
                    
                    if crosses_date_line:
                        # For crossing date line: longitude should be >= min_lon OR <= max_lon
                        lon_in_range = lon >= min_lon or lon <= max_lon
                    else:
                        # Normal case: longitude should be between min_lon and max_lon
                        lon_in_range = min_lon <= lon <= max_lon
                    
                    if lat_in_range and lon_in_range:
                        sites_within_bbox[site_name] = {
                            'lat': lat,
                            'lon': lon,
                            'elevation': elevation
                        }
                        
                except (ValueError, IndexError) as e:
                    print(f"Error parsing line: {line}")
                    print(f"Error details: {e}")
                    continue
            
            #print(f"Found {len(sites_within_bbox)} AERONET sites within the specified bounding box")
            return sites_within_bbox
            
        except Exception as e:
            print(f"Error downloading or parsing AERONET site list: {e}")
            return {}
    

    
    def download_aeronet_data(self, site_name, variables):
        """
        Download AERONET data for a specific site and date range
        
        Parameters:
        -----------
        site_name : str
            AERONET site name
        
        variables: list
            Variable names to download (e.g., ['AOD_500nm'] or ['AOD_500nm', 'AOD_380nm']
        
        Returns:
        --------
        str or None
            Path to downloaded CSV file, or None if download failed
        """
        print('     ')
        print(f"Searching AERONET data for {site_name}...")
        
        # Convert dates to the format needed for the AERONET website
        beg_year = self.start_date.year
        beg_month = self.start_date.month
        beg_day = self.start_date.day
        
        end_year = self.end_date.year
        end_month = self.end_date.month
        end_day = self.end_date.day
        
        
        # The base URL of AERONET data for the Version 3 data
        base_url = "https://aeronet.gsfc.nasa.gov/cgi-bin/print_web_data_v3"
        
        # Prepare payload
        payload = {
            "site": site_name,
            "year": beg_year,
            "month": beg_month,
            "day": beg_day,
            "year2": end_year,
            "month2": end_month,
            "day2": end_day,
            self.aeronet_data_type: 1,
            "AVG": self.aeronet_data_format
        }
        
        try:
            # Web scraping the data
            response = requests.get(base_url, params=payload, timeout=30)
            
            print(f" Loop through--> {response.url}")
            
            if response.status_code != 200:
                print(f"The url:")
                print(f"--> {response.url}")
                print(f"Is not reachable. Please check your settings.")
                return None
                       
            # Find data start
            lines = response.text.split('\n')
            header_idx = None
            for i, line in enumerate(lines):
                if 'Date(dd:mm:yyyy)' in line:
                    header_idx = i
                    break
            
            # FIX: Convert decimal dates to proper format
            from datetime import datetime, timedelta
            for i in range(header_idx + 1, len(lines)):
                if lines[i].strip() and ',' in lines[i]:
                    parts = lines[i].split(',')
                    try:
                        # If first column is decimal, convert it
                        if '.' in parts[0] and float(parts[0]) > 0:
                            excel_date = float(parts[0])
                            converted_date = datetime(1899, 12, 30) + timedelta(days=excel_date)
                            parts[0] = converted_date.strftime('%d:%m:%Y')
                            lines[i] = ','.join(parts)
                    except:
                        pass  # Keep original if conversion fails
            
            # Create DataFrame
            data_text = '\n'.join(lines[header_idx:])
            df = pd.read_csv(io.StringIO(data_text))
            df.columns = df.columns.str.strip()
            
            # Keep date columns + requested variables  
            date_cols = [col for col in df.columns if 'date' in col.lower() or 'time' in col.lower()]
            keep_cols = date_cols + [var for var in variables if var in df.columns]
            
            df_filtered = df[keep_cols]
            
            # Save to CSV
            # Create output filename for this site's data
            data_path = os.path.join(self.output_dir, 'data/AERONET')
            os.makedirs(data_path, exist_ok=True)
            tmp_csv_file = f'{data_path}/AERONET.{site_name}.{beg_year}{beg_month:02d}{beg_day:02d}.{end_year}{end_month:02d}{end_day:02d}.AVG{self.aeronet_data_format}.{self.aeronet_data_type}.csv'
            
            df_filtered.to_csv(tmp_csv_file, index=False)
            
            # Check if file has valid data
            print('os.path.getsize(tmp_csv_file)', os.path.getsize(tmp_csv_file))
            if os.path.getsize(tmp_csv_file) > 100:  # Using smaller threshold since we're filtering to just one site
                print(f"Downloaded AERONET data for {site_name}")
                return tmp_csv_file
            else:
                print(f"Sorry, no data available for {site_name} in the specified period")
                os.remove(tmp_csv_file)
                return None
                
        except Exception as e:
            print(f"Error downloading AERONET data for {site_name}: {str(e)}")
            return None
        
    
    def read_aeronet_data(self, filename):
        """
        Read and parse AERONET CSV data and find AOD column closest to target wavelength
        
        Parameters:
        -----------
        filename : str
            Path to AERONET CSV file
        target_wavelength : int
            Target wavelength for AOD data (default: 500 nm)
            
        Returns:
        --------
        tuple
            (DataFrame with datetime and aod columns, actual_wavelength used)
            Returns (None, None) if reading fails
        """
        
        try:
            # Since we now write only header + data rows, we can read directly
            # without needing to search for the header line
            df = pd.read_csv(filename)
            
            # Check if the expected columns exist
            if 'Date(dd:mm:yyyy)' not in df.columns or 'Time(hh:mm:ss)' not in df.columns:
                print("Required date/time columns not found in AERONET file")
                print(f"Available columns: {df.columns.tolist()}")
                return None, None
            
            # Convert date/time columns to datetime
            df['datetime'] = pd.to_datetime(
                df['Date(dd:mm:yyyy)'].astype(str) + ' ' + df['Time(hh:mm:ss)'].astype(str),
                format='%d:%m:%Y %H:%M:%S',
                errors='coerce'
            )
            
            print( df['datetime'])
            # Find all AOD columns and extract wavelengths
            aod_columns = [col for col in df.columns if 'AOD' in col and any(char.isdigit() for char in col)]
            
            if not aod_columns:
                print("No AOD columns found in AERONET data")
                print(f"Available columns: {df.columns.tolist()}")
                return None, None
            
            #print(f"Found AOD columns: {aod_columns}")
            
            # Extract wavelengths from column names
            wavelengths = []
            for col in aod_columns:
                # Extract numbers from column name
                numbers = [int(s) for s in col.split('_') if s.isdigit()]
                if numbers:
                    wavelengths.append(numbers[0])
                else:
                    # Try to extract from other formats like "AOD_500nm"
                    import re
                    match = re.search(r'(\d+)', col)
                    if match:
                        wavelengths.append(int(match.group(1)))
            
            if not wavelengths:
                print("Could not extract wavelengths from AOD column names")
                return None, None
            
            # Find the wavelength closest to target
            wavelength_diffs = [abs(w - self.target_wavelength) for w in wavelengths]
            closest_idx = wavelength_diffs.index(min(wavelength_diffs))
            selected_wavelength = wavelengths[closest_idx]
            selected_aod_col = aod_columns[closest_idx]
            
            print(f"Selected AOD column: {selected_aod_col} (wavelength: {selected_wavelength}nm)")
            
            # Create working dataframe with selected AOD column
            df['aod'] = pd.to_numeric(df[selected_aod_col], errors='coerce')
            
            # Remove rows with NaN datetime or AOD values
            df = df.dropna(subset=['aod', 'datetime'])
            
            # Remove invalid AOD values
            df = df[(df['aod'] >= 0) & (df['aod'] <= 5.0)]
            
            print(f"Loaded {len(df)} AERONET measurements using {selected_aod_col}")
            
            # Return both the dataframe and the actual wavelength used
            return df[['datetime', 'aod']].copy(), selected_wavelength
                
        except Exception as e:
            print(f"Error reading AERONET file {filename}: {str(e)}")
            
            # Debug information - show first few lines of the file
            try:
                with open(filename, 'r') as f:
                    print("First 10 lines of the file:")
                    for i, line in enumerate(f.readlines()[:10]):
                        print(f"Line {i}: {line.strip()}")
            except:
                pass
            
            return None, None
        
    def download_tropomi_data(self):
        """
        Search and download TROPOMI data within a specified bounding box
        
            
        Returns:
        --------
        data_dir: str
            the directory where the data are downloaded
        """
        print(f"Searching and downloading TROPOMI data within the specified bounding box...")
        
        data_path = None  # Initialize data_path
        
        try:
        
            # First time setup - this will prompt for credentials and save them
            #earthaccess.login(persist=True)

            # After first setup, subsequent logins will use saved credentials
            auth = earthaccess.login(strategy="interactive", persist=True)
            
            min_lon = self.bounding_box['min_lon'] 
            max_lon = self.bounding_box['max_lon']
            min_lat = self.bounding_box['min_lat']
            max_lat = self.bounding_box['max_lat']
            
            # To download multiple files, change the second temporal parameter
            results = earthaccess.search_data(
                short_name=self.tropomi_data_shortname,
                version=self.tropomi_data_version,
                temporal=(self.start_date, self.end_date), 
                bounding_box=(min_lon, min_lat, max_lon, max_lat)
            )
            for result in results:
                data_urls = result['umm']['RelatedUrls']
                # Filter for data URLs (not browse/metadata URLs)
                data_links = [url['URL'] for url in data_urls if url['Type'] == 'GET DATA']
                print("  Data URLs:", data_links)
    
            print(f"  Downloading... It may take longer for more files")
            # Download granules to local path
            #data_path = os.path.join(self.output_dir, 'data/'+self.tropomi_data_shortname+'.'+self.tropomi_data_version) 
            data_path = os.path.join(self.output_dir, 'data', f"{self.tropomi_data_shortname}.{self.tropomi_data_version}")
            downloaded_files = earthaccess.download(
                results,
                local_path=data_path, # Change this string to download to a different path
                )
            return data_path
         
        except Exception as e:
            print(f"Error searching or downloading TROPOMI data: {str(e)}")
            print(f"Error type: {type(e).__name__}")
            if data_path:
                print(f"Attempted data_path: {data_path}")
            import traceback
            traceback.print_exc()
            return None
            
    def read_tropomi_data(self, filename):
        """
        Read TROPOMI aerosol data from netCDF file
        
        Parameters:
        -----------
        filename : str
            Path to TROPOMI netCDF file
        target_wavelength : int
            Target wavelength for AOD data (default: 500 nm)
        """      
        try:
            with nc.Dataset(filename, 'r') as ds:
                # Read geolocation data from GEODATA group
                lat = ds.groups['GEODATA']['latitude'][:]
                lon = ds.groups['GEODATA']['longitude'][:]
                
                # Read 1D time data in original data 
                delta_time = ds.groups['GEODATA']['delta_time'][:]
                
                
                # Parse reference time from delta_time units
                time_units = ds.groups['GEODATA']['delta_time'].units
                # Extract reference time from units string: "milliseconds since 2025-01-08 00:00:00"
                ref_time_str = time_units.split('since ')[1]
                ref_time = datetime.strptime(ref_time_str, '%Y-%m-%d %H:%M:%S')
                      
                # Calculate satellite overpass time for the full scene (will be refined in spatial_collocation)
                # Convert delta_time from milliseconds to seconds and add to reference time            
                mean_delta_seconds = np.nanmean(delta_time) / 1000.0
                max_delta_seconds = np.nanmax(delta_time) / 1000.0
                min_delta_seconds = np.nanmin(delta_time) / 1000.0
                
                satellite_time_min = ref_time + timedelta(seconds=min_delta_seconds)
                satellite_time_max = ref_time + timedelta(seconds=max_delta_seconds)
                satellite_time_mean = ref_time + timedelta(seconds=mean_delta_seconds)
                         
                
                #print(f"satellite_time_min: {satellite_time_min}")
                #print(f"satellite_time_max: {satellite_time_max}")
                #print(f"satellite_time_mean: {satellite_time_mean}")
        
                
                # Mean of satellite overpass time for all pixels
                #satellite_time = ref_time + timedelta(seconds=mean_delta_seconds)
                satellite_time = satellite_time_mean
                
                # Read aerosol optical depth from SCIDATA group
                # FinalAerosolOpticalDepth has dimensions (scanline, ground_pixel, Wavelengths)
                aod_all_wavelengths = ds.groups['SCIDATA']['FinalAerosolOpticalDepth'][:]
                
                # Get wavelength information
                wavelengths = ds['Wavelengths'][:]
                print(f"Available wavelengths: {wavelengths} nm")
                
                # Select wavelength closest to yarget wavelength, e.g., 380nm (typically used for validation)
                #target_wavelength = 380.0
                wl_diff = np.abs(wavelengths - self.target_wavelength)
                best_wl_idx = np.argmin(wl_diff)
                selected_wavelength = wavelengths[best_wl_idx]
                
                print(f"Selected wavelength: {selected_wavelength} nm (closest to {self.target_wavelength} nm)")
                
                # Extract AOD for selected wavelength
                aod = aod_all_wavelengths[:, :, best_wl_idx]
                
                # Read quality flags (FinalAlgorithmFlags)
                qa_flags = ds.groups['SCIDATA']['FinalAlgorithmFlags'][:]
                
                # Create quality mask
                # Valid AOD values: not fill value, within reasonable range, and finite
                aod_fill_value = ds.groups['SCIDATA']['FinalAerosolOpticalDepth']._FillValue
                lat_fill_value = ds.groups['GEODATA']['latitude']._FillValue
                lon_fill_value = ds.groups['GEODATA']['longitude']._FillValue
                
                # Create comprehensive quality mask
                valid_mask = (
                    (aod != aod_fill_value) & 
                    (lat != lat_fill_value) & 
                    (lon != lon_fill_value) &
                    (aod > 0) & 
                    (aod < 5.0) & 
                    np.isfinite(aod) & 
                    np.isfinite(lat) & 
                    np.isfinite(lon) &
                    (lat >= -90) & (lat <= 90) &
                    (lon >= -180) & (lon <= 180)
                )
                
                # Additional quality filtering based on FinalAlgorithmFlags if needed
                # For now, we'll accept all flag values (0) as they indicate data source
                qa_fill_value = ds.groups['SCIDATA']['FinalAlgorithmFlags']._FillValue
                valid_qa_mask = (qa_flags != qa_fill_value) & (qa_flags == 0)
    
                # Additional quality filtering based on AIRSCO_Flags if needed
                # For now, we'll accept all flag values (1, 2, 3) as they indicate data source
                # qa_fill_value = ds.groups['SCIDATA']['AIRSCO_Flags']._FillValue
                # valid_qa_mask = (qa_flags != qa_fill_value) & (qa_flags >= 1) & (qa_flags <= 3)
                
                # Combine quality masks
                final_mask = valid_mask & valid_qa_mask
                
                # Store original shape before flattening
                original_shape = lat.shape  # This should be (4173, 450)
                
                # Flatten arrays and apply mask
                lat_flat = lat[final_mask]
                lon_flat = lon[final_mask]
                aod_flat = aod[final_mask]
                qa_flat = qa_flags[final_mask]
                
                
                #print(f"Loaded TROPOMI data: {len(aod_flat)} valid pixels out of {lat.size} total pixels")
                #print(f"Data quality distribution: {np.unique(qa_flat, return_counts=True)}")
                #print(f"Satellite overpass time (mean of all pixles): {satellite_time}")
                #print(f"AOD range: {np.min(aod_flat):.3f} - {np.max(aod_flat):.3f}")
                
                # Convert delta_time to 2D by broadcasting across ground pixels
                # Each scanline has the same time for all ground pixels
                delta_time_2d = np.broadcast_to(delta_time[:, np.newaxis], lat.shape)
            
                # print(f"DEBUG - delta_time shape: {delta_time.shape}")
                # print(f"DEBUG - delta_time_2d shape: {delta_time_2d.shape}")
                # print(f"DEBUG - lat shape: {lat.shape}")
                
                delta_time_flat = delta_time_2d[final_mask]  # Now flattened consistently
                
                return {
                    'latitude': lat_flat,
                    'longitude': lon_flat,
                    'aod': aod_flat,
                    'qa_flags': qa_flat,
                    'wavelength': selected_wavelength,
                    'satellite_time': satellite_time,
                    #'delta_time': delta_time,  # Include the original delta_time array
                    'delta_time': delta_time_flat,  # Now a flattened array, same length as lat/lon
                    'ref_time': ref_time,      # Include the reference time
                    'original_shape': original_shape,  # Add this line
                    'filename': os.path.basename(filename)
                }
                
        except Exception as e:
            print(f"Error reading TROPOMI file {filename}: {str(e)}")
            print("Attempting to diagnose file structure...")
            
            # Diagnostic information
            try:
                with nc.Dataset(filename, 'r') as ds:
                    print(f"Root groups: {list(ds.groups.keys())}")
                    if 'SCIDATA' in ds.groups:
                        print(f"SCIDATA variables: {list(ds.groups['SCIDATA'].variables.keys())}")
                    if 'GEODATA' in ds.groups:
                        print(f"GEODATA variables: {list(ds.groups['GEODATA'].variables.keys())}")
                    print(f"Root variables: {list(ds.variables.keys())}")
            except:
                pass
            
            return None
    
    def spatial_collocation(self, tropomi_data, site_coords):
        """
        Find TROPOMI pixels within spatial radius of AERONET site
        """
        # Calculate distances
        
        print('site_coords:', site_coords)
        distances = self.haversine_distance(
            tropomi_data['latitude'], tropomi_data['longitude'],
            site_coords['lat'], site_coords['lon']
        )
        
        # Find pixels within radius
        within_radius = distances <= self.spatial_radius_km
        
        # Print summary information about within_radius
        total_pixels = len(within_radius)
        matched_pixels = np.sum(within_radius)
        percent_matched = (matched_pixels / total_pixels) * 100 if total_pixels > 0 else 0

        print(f"\nSpatial collocation results:")
        print(f"Found {total_pixels} TROPOMI pixels within {self.spatial_radius_km} km")
        print(f"  Total TROPOMI pixels analyzed: {total_pixels:,}")
        print(f"  Pixels within {self.spatial_radius_km} km radius: {matched_pixels:,} ({percent_matched:.2f}%)")
        
        if np.any(within_radius):
            # Get the original dimensions from the data structure
            # TROPOMI data: latitude(scanline=4173, ground_pixel=450), delta_time(scanline=4173, ground_pixel=450). 
            #               Where the scanline serves as the index for measurements in the along-track direction; and each of these across-track measurements is a ground pixel.
            # The flattened arrays need to be mapped back to scanline indices
            
            #original_shape = tropomi_data.get('original_shape', (4173, 450))  # Default TROPOMI dimensions
            n_scanlines, n_ground_pixels = tropomi_data["original_shape"]
            
            # Get indices of collocated pixels in the flattened array
            flat_indices = np.where(within_radius)[0]
            
            # Print first few and last few indices
            max_display = 10  # Number of elements to display from start/end

            if len(flat_indices) <= max_display * 2:
                # If array is small, print all indices
                indices_str = ', '.join(map(str, flat_indices))
            else:
                # Print first and last few indices
                first_indices = ', '.join(map(str, flat_indices[:max_display]))
                last_indices = ', '.join(map(str, flat_indices[-max_display:]))
                indices_str = f"{first_indices}, ..., {last_indices}"
         
            # Get delta_time values for the corresponding scanlines
            collocated_delta_time = tropomi_data['delta_time'][flat_indices]
            
            # Calculate satellite overpass time using only the collocated pixels
            # Convert delta_time from milliseconds to seconds and add to reference time            
            mean_delta_seconds = np.nanmean(collocated_delta_time) / 1000.0
            max_delta_seconds = np.nanmax(collocated_delta_time) / 1000.0
            min_delta_seconds = np.nanmin(collocated_delta_time) / 1000.0
            
            satellite_time_min = tropomi_data['ref_time'] + timedelta(seconds=min_delta_seconds)
            satellite_time_max = tropomi_data['ref_time'] + timedelta(seconds=max_delta_seconds)
            satellite_time_mean = tropomi_data['ref_time'] + timedelta(seconds=mean_delta_seconds)
            
            for i in range(len(flat_indices)):
                # Convert delta_time (seconds) to timedelta before adding to ref_time
                delta_seconds = tropomi_data['delta_time'][flat_indices[i]] / 1000.0  # Convert from milliseconds to seconds
                satellite_time = tropomi_data['ref_time'] + pd.Timedelta(seconds=delta_seconds)
                
                # Get actual coordinate values
                pixel_lat = tropomi_data['latitude'][flat_indices[i]]
                pixel_lon = tropomi_data['longitude'][flat_indices[i]]
                pixel_distance = distances[flat_indices[i]]
                
                # Get variable value
                aod_values = tropomi_data['aod'][flat_indices[i]]

                print(f"    Flat index {flat_indices[i]}: lat={pixel_lat:8.4f}°, lon={pixel_lon:9.4f}°, "
                      f"dist={pixel_distance:6.2f}km, {satellite_time.strftime('%H:%M:%S.%f')[:-3]}, "
                      f"value={aod_values:6.2f} ")
            
            # Calculate the AOD statistics for matched pixels
            aod_values = tropomi_data['aod'][flat_indices]
            print(f"\n  AOD statistics for matched pixels:")
            print(f"    Mean: {np.mean(aod_values):.4f}")
            print(f"    Median: {np.median(aod_values):.4f}")
            print(f"    Min: {np.min(aod_values):.4f}")
            print(f"    Max: {np.max(aod_values):.4f}")
            print(f"    Std dev: {np.std(aod_values):.4f}")
            
            # Recalculate satellite time based on collocated pixels
            #satellite_time = tropomi_data['ref_time'] + timedelta(seconds=mean_delta_seconds)
            satellite_time = satellite_time_mean
            
            collocated_data = {
                'latitude': tropomi_data['latitude'][within_radius],
                'longitude': tropomi_data['longitude'][within_radius],
                'aod': tropomi_data['aod'][within_radius],
                'distance_km': distances[within_radius],
                'satellite_time': satellite_time,
                'wavelength': tropomi_data['wavelength'],
                'filename': tropomi_data['filename']
            }
            
            #print(f"Found {len(collocated_data['aod'])} TROPOMI pixels within {self.spatial_radius_km} km")
            print(f"collocated pixels: satellite_time_min: {satellite_time_min}")
            print(f"collocated pixels: satellite_time_max: {satellite_time_max}")
            print(f"collocated pixels: satellite_time_mean: {satellite_time_mean}")
            print(f"Satellite overpass time (mean of collocated pixels): {satellite_time}")
            return collocated_data
        else:
            return None
        
    
    def temporal_collocation(self, tropomi_data, aeronet_data):
        """
        Match satellite with the averaged AERONET measurements within temporal window of satellite overpass time 
        """
        satellite_time = tropomi_data['satellite_time']
        
        # Calculate time differences in hours
        time_diff = np.abs((aeronet_data['datetime'] - satellite_time).dt.total_seconds() / 3600)
        
        # Find measurements within temporal window
        within_window = time_diff <= self.temporal_window_hours
        
        # Print detailed information about temporal matching
        total_measurements = len(time_diff)
        measurements_in_window = np.sum(within_window)

        print(f"\nTemporal collocation analysis:")
        print(f"  Satellite overpass time: {satellite_time}")
        print(f"  Temporal window: ±{self.temporal_window_hours} hours")
        print(f"  Total AERONET measurements available: {total_measurements}")
        print(f"  Measurements within temporal window: {measurements_in_window}")
        
        if np.any(within_window):
            
            # Get all AERONET measurements within the temporal window
            windowed_aeronet = aeronet_data.loc[within_window].copy()
            
            # Show all data for each measurement                    
            # Show the actual data
            print(f"\nData preview of windowed_aeronet:")
            print(windowed_aeronet)
            
            # Calculate average of all measurements within window            
            # Create DataFrame with explicit single row
            averaged_aeronet = pd.DataFrame({
                'datetime': [satellite_time],  # Use list to create single row
                'aeronet_aod_mean': [windowed_aeronet['aod'].mean()],
                'aeronet_aod_std': [windowed_aeronet['aod'].std()],
                'n_measurements': [np.sum(within_window)],
                'time_window_start': [windowed_aeronet['datetime'].min()],
                'time_window_end': [windowed_aeronet['datetime'].max()]
            })
            
            print(f"\n averaged_aeronet:")
            print(averaged_aeronet)

            print(f"Found {np.sum(within_window)} AERONET measurements within ±{self.temporal_window_hours} hours")
            print(f"Averaged measurements with mean time difference: {time_diff[within_window].mean():.3f} ± {time_diff[within_window].std():.3f} hours")

            return averaged_aeronet
        else:
            return None
    
    def collocate_data(self):
        """
        Main collocation function that processes TROPOMI files and matches with AERONET data
        
        Parameters:
        -----------
        tropomi_files : list
            List of TROPOMI file paths
        start_date : datetime
            Start date for analysis
        end_date : datetime  
            End date for analysis
        target_wavelength : int
            Target wavelength for AERONET AOD comparison
        """
        
        all_collocations = []
        
        # Download AERONET data for all sites
        aeronet_data = {}
        for site_name, coords in self.aeronet_sites.items():
            aeronet_file = self.download_aeronet_data(site_name, ['AOD_'+str(self.target_wavelength)+'nm'])
            if aeronet_file:
                aeronet_result = self.read_aeronet_data(aeronet_file)
                if aeronet_result[0] is not None:
                    aeronet_data[site_name] = aeronet_result[0]  # Store just the DataFrame
                    
                    actual_wavelength = aeronet_result[1]        # Get the wavelength
                    print(f"Using AERONET AOD at {actual_wavelength}nm for {site_name}")
                else:
                    continue
        
        print(f"\nSuccessfully loaded AERONET data for {len(aeronet_data)} sites")
        
        # Process each TROPOMI file
        
        for tropomi_file in self.tropomi_files:
            print(f"\nProcessing TROPOMI file: {os.path.basename(tropomi_file)}")
            
            tropomi_data = self.read_tropomi_data(tropomi_file)
            if tropomi_data is None:
                continue
            
            # Check if TROPOMI data is within our analysis period
            if not (self.start_date <= tropomi_data['satellite_time'] <= self.end_date):
                print(f"TROPOMI data outside analysis period, skipping...")
                continue
            
            # Collocate with each AERONET site
            for site_name, site_coords in self.aeronet_sites.items():
                if site_name not in aeronet_data:
                    continue
                
                print(f"\nCollocating with {site_name}...")
                
                # The 1st step: Spatial collocation first then send the result to do the temporal collocation
                spatial_match = self.spatial_collocation(tropomi_data, site_coords)
                print(f" Spatial collocation is done")
                if spatial_match is None:
                    continue
                
                # The 2nd step: Temporal collocation with target wavelength
                temporal_match = self.temporal_collocation(spatial_match, aeronet_data[site_name])
                print(f" Temporal collocation is done")
                if temporal_match is None:
                    continue
                
                # Store collocation results
                for _, aeronet_row in temporal_match.iterrows():
                    # Average TROPOMI AOD within spatial radius
                    tropomi_aod_mean = np.mean(spatial_match['aod'])
                    tropomi_aod_std = np.std(spatial_match['aod'])
                    
                    # Add this before the problematic line to see available columns
                    print("Available AERONET columns:", aeronet_row.index.tolist())
                    print(f'Look for aeronet_aod_{actual_wavelength}nm' ) 
                    
                    collocation= {
                        'site_name': site_name,
                        'site_lat': site_coords['lat'],
                        'site_lon': site_coords['lon'],
                        'aeronet_datetime': temporal_match['datetime'].iloc[0],  # This is now the collocated time
                        'aeronet_aod_mean': temporal_match['aeronet_aod_mean'].iloc[0],  
                        'aeronet_aod_std': temporal_match['aeronet_aod_std'].iloc[0] if 'aeronet_aod_std' in temporal_match.columns else None,
                        'aeronet_measurement_count': temporal_match['n_measurements'].iloc[0],
                        'time_window_start': temporal_match['time_window_start'].iloc[0],
                        'time_window_end': temporal_match['time_window_end'].iloc[0],
                        'tropomi_datetime': spatial_match['satellite_time'],                   
                        'tropomi_aod_mean': tropomi_aod_mean,
                        'tropomi_aod_std': tropomi_aod_std,
                        'tropomi_pixel_count': len(spatial_match['aod']),
                        'tropomi_wavelength': spatial_match['wavelength'],  # Fixed: Now available in spatial_match
                        'min_distance_km': np.min(spatial_match['distance_km']),
                        'max_distance_km': np.max(spatial_match['distance_km']),
                        'tropomi_filename': spatial_match['filename'],

                    } 
        
                    all_collocations.append(collocation)
            
        # Convert to DataFrame
        if all_collocations:
            df_collocations = pd.DataFrame(all_collocations)
            # Sort the DataFrame by site_name first, then by aeronet_datetime
            df_collocations = df_collocations.sort_values(by=['site_name', 'aeronet_datetime'])
    
            print(f"\nTotal collocations found: {len(df_collocations)}")
            
            # Save results
            output_file = f"{self.output_dir}/data/collocation_results_{self.target_wavelength}nm.csv"
            df_collocations.to_csv(output_file, index=False)
            print(f"Results saved to: {output_file}")
            
            return df_collocations
        else:
            print("\nNo collocations found!")
            return None
    
    def plot_collocation_map(self, df_collocations):
        """
        Plot collocation locations on a map
        """
        if df_collocations is None or len(df_collocations) == 0:
            print("No data to plot")
            return
        
        fig = plt.figure(figsize=(12, 8))
        ax = plt.axes(projection=ccrs.PlateCarree())
        
        # Set extent for the selected area
        min_lon = self.bounding_box['min_lon'] 
        max_lon = self.bounding_box['max_lon']
        min_lat = self.bounding_box['min_lat']
        max_lat = self.bounding_box['max_lat']
        
        ax.set_extent([min_lon, max_lon, min_lat, max_lat], ccrs.PlateCarree())
        
        # Add map features
        ax.add_feature(cfeature.COASTLINE)
        ax.add_feature(cfeature.BORDERS)
        ax.add_feature(cfeature.STATES)
        ax.add_feature(cfeature.OCEAN, alpha=0.3)
        ax.add_feature(cfeature.LAND, alpha=0.3)
        
        # Plot AERONET sites
        #for site_name, coords in self.la_aeronet_sites.items():
        #    ax.plot(coords['lon'], coords['lat'], 'rs', markersize=8, 
        #           transform=ccrs.PlateCarree(), label='AERONET Sites' if site_name == list(self.la_aeronet_sites.keys())[0] else "")
        #    ax.text(coords['lon'], coords['lat'] + 0.05, site_name, 
        #           transform=ccrs.PlateCarree(), ha='center', fontsize=8)
        
        # Plot collocation points
        unique_sites = df_collocations.groupby('site_name').first()
        for site_name, row in unique_sites.iterrows():
            site_collocations = df_collocations[df_collocations['site_name'] == site_name]
            ax.plot(row['site_lon'], row['site_lat'], 'bo', markersize=6, 
                   transform=ccrs.PlateCarree(), alpha=0.7)
            # Site name
            ax.text(row['site_lon'], row['site_lat'] + 0.05, site_name,
                    transform=ccrs.PlateCarree(), ha='center', fontsize=8)
            
            # Draw circles showing spatial radius
            circle = plt.Circle((row['site_lon'], row['site_lat']), 
                              self.spatial_radius_km/111.32, # Convert km to degrees (approximate)
                              fill=False, linestyle='--', alpha=0.5, 
                              transform=ccrs.PlateCarree())
            ax.add_patch(circle)
        
        ax.gridlines(draw_labels=True, alpha=0.5)
        ax.set_title(f'TROPOMI-AERONET Collocations in the Area of Interest\n'
  #                  f'Total Number of Collocated Observations: {len(df_collocations)}', fontsize=14)
                   f'Total Number of Collocated Sites: {len(unique_sites)}', fontsize=14)
        
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"{self.output_dir}/plots/collocation_map_{int(self.spatial_radius_km)}km.png", dpi=300, bbox_inches='tight')
        plt.show()
    

    def plot_site_time_series(self, df_results, site_name, tropomi_wavelength, save_plots=True):
        """
        Plot time series of collocated TROPOMI and AERONET AOD data for a specific site
        
        Parameters:
        -----------
        df_results : pandas.DataFrame
            DataFrame containing collocation results
        site_name : str
            Name of the AERONET site to plot
        target_wavelength : int
            Target wavelength for AERONET data
        tropomi_wavelength : float
            Actual wavelength used from TROPOMI data in collocation(from read_tropomi_data)
        save_plots : bool
            Whether to save plots to file
        """
        # Filter data for the specific site
        site_data = df_results[df_results['site_name'] == site_name].copy()
        
        if len(site_data) == 0:
            print(f"No data found for site: {site_name}")
            return None, None  # Return consistent tuple format
        
        # Convert datetime strings to datetime objects if needed
        if isinstance(site_data['aeronet_datetime'].iloc[0], str):
            site_data['aeronet_datetime'] = pd.to_datetime(site_data['aeronet_datetime'])
        
        # Sort by time
        site_data = site_data.sort_values('aeronet_datetime')
        
        # Create the plot
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Get site location information - FIXED VERSION
        site_location = ""
        try:
            if hasattr(self, 'aeronet_sites') and site_name in self.aeronet_sites:
                site_coords = self.aeronet_sites[site_name]
                
                # Handle dictionary format
                if isinstance(site_coords, dict):
                    lat = site_coords.get('lat')
                    lon = site_coords.get('lon')
                    elevation = site_coords.get('elevation')
                    
                    if lat is not None and lon is not None:
                        # Format: (lat: 34.14°, lon: -118.13°, Alt: 260.0m)
                        site_location = f" (lat: {lat:.2f}°, lon: {lon:.2f}°"
                        
                        # Add elevation if available
                        if elevation is not None:
                            site_location += f", Alt: {elevation:.1f}m)"
                        else:
                            site_location += ")"
                    else:
                        print(f"Warning: Could not find lat/lon in dictionary for {site_name}")
                        
                elif isinstance(site_coords, (list, tuple)):
                    if len(site_coords) >= 2:
                        lat, lon = site_coords[0], site_coords[1]
                        site_location = f" (lat: {lat:.2f}°, lon: {lon:.2f}°)"
                    else:
                        print(f"Warning: Insufficient coordinate data for {site_name}")
                else:
                    print(f"Warning: Unexpected coordinate format for {site_name}: {type(site_coords)}")
                    
        except Exception as e:
            print(f"Warning: Could not extract coordinates for {site_name}: {str(e)}")
            # Try to get coordinates from the DataFrame instead
            if len(site_data) > 0 and 'site_lat' in site_data.columns and 'site_lon' in site_data.columns:
                lat = site_data['site_lat'].iloc[0]
                lon = site_data['site_lon'].iloc[0]
                site_location = f" (lat: {lat:.2f}°, lon: {lon:.2f}°)"
        
        # Get AERONET data information from the DataFrame if available
        aeronet_data_info = ""
              
        aeronet_data_info = f" - {self.aeronet_data_type} AVG{self.aeronet_data_format}".strip(" -")
        #print(f"aeronet_data_info is {aeronet_data_info}")
        
        
        # Plot AERONET AOD data (all available points) - use target_wavelength
        ax.scatter(site_data['aeronet_datetime'], site_data['aeronet_aod_mean'], 
                  color='black', marker='o', s=50, alpha=0.7, 
                  label=f'AERONET AOD ({self.target_wavelength}nm - {aeronet_data_info})', zorder=3)
        
        # Plot TROPOMI AOD data (collocated points) - use tropomi_wavelength
        ax.scatter(site_data['aeronet_datetime'], site_data['tropomi_aod_mean'], 
                  color='dodgerblue', marker='s', s=60, alpha=0.8, 
                  #label=f'TROPOMI AOD ({tropomi_wavelength:.1f}nm - collocated)', zorder=3)
                  label=f'TROPOMI AOD ({int(tropomi_wavelength)}nm - collocated)', zorder=3)
        
        # Add error bars for TROPOMI if standard deviation is available
        if 'tropomi_aod_std' in site_data.columns:
            ax.errorbar(site_data['aeronet_datetime'], site_data['tropomi_aod_mean'],
                       yerr=site_data['tropomi_aod_std'], fmt='none', 
                       color='blue', alpha=0.5, capsize=3, zorder=2, 
                       elinewidth=0.1,  # Makes the error bar lines thinner
                       capthick=0.5,    # Makes the caps thinner too
                       )
            
        # Add error bars for AERONET if standard deviation is available
        if 'aeronet_aod_std' in site_data.columns:
            ax.errorbar(site_data['aeronet_datetime'], site_data['aeronet_aod_mean'],
                       yerr=site_data['aeronet_aod_std'], fmt='none', 
                       color='gray', alpha=0.5, capsize=3, zorder=2, 
                       elinewidth=0.1,  # Makes the error bar lines thinner
                       capthick=0.5,    # Makes the caps thinner too
                       )
        
        # Connect points with lines for better visualization
        ax.plot(site_data['aeronet_datetime'], site_data['aeronet_aod_mean'], 
               color='black', alpha=0.3, linewidth=1, zorder=1)
        ax.plot(site_data['aeronet_datetime'], site_data['tropomi_aod_mean'], 
               color='dodgerblue', alpha=0.3, linewidth=1, zorder=1)
        
        # Formatting
        ax.set_xlabel('Date and time (UTC)', fontsize=12)
        ax.set_ylabel('Aerosol Optical Depth', fontsize=12)
        
        # Get time range
        start_date = site_data['aeronet_datetime'].min()
        end_date = site_data['aeronet_datetime'].max()
        time_range = f"{start_date.strftime('%Y-%m-%dT%H')} - {end_date.strftime('%Y-%m-%dT%H')}"
  
        # Get matching criteria information
        matching_info = f"{self.spatial_radius_km}km and ±{self.temporal_window_hours}hr"
        
        # Updated title with site location
        ax.set_title(f'AERONET AOD Time Series - Matching TROPOMI within {matching_info} \nSite: {site_name}{site_location} \nTime_range:{time_range}', 
                    fontsize=14, fontweight='bold')
        
        # Position legend to avoid overlap with stats text (top right instead of default)
        ax.legend(fontsize=11, loc='upper right')
        ax.grid(True, alpha=0.3)
        
        # Format x-axis dates to show year but not seconds
        from matplotlib.dates import DateFormatter
        #date_formatter = DateFormatter('%Y-%m-%d')  # Shows year-month-day, no time
        date_formatter = DateFormatter('%Y-%m-%dT%H')  # 2023-08-15 14:30
        ax.xaxis.set_major_formatter(date_formatter)
        # More control over label rotation
        plt.xticks(rotation=45, ha='right')  # Rotate 45 degrees, align to right

        fig.autofmt_xdate()
        
        # Add statistics text box (positioned at top left to avoid legend)
        correlation = np.corrcoef(site_data['aeronet_aod_mean'], site_data['tropomi_aod_mean'])[0,1]
        rmse = np.sqrt(np.mean((site_data['aeronet_aod_mean'] - site_data['tropomi_aod_mean'])**2))
        bias = np.mean(site_data['tropomi_aod_mean'] - site_data['aeronet_aod_mean'])
        n_points = len(site_data)
        
        stats_text = f'N = {n_points}\nR = {correlation:.3f}\nRMSE = {rmse:.3f}\nBias = {bias:.3f}'
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        plt.tight_layout()
        
        # Save plot if requested
        if save_plots:
            plot_filename = f"time_series_{site_name.replace(' ', '_')}.png"
            plot_path = os.path.join(self.output_dir, 'plots/'+plot_filename)
            plt.savefig(plot_path, dpi=300, bbox_inches='tight')
            print(f"     Time series plot saved for {site_name}: {plot_path}")
        
        plt.show()
        return fig, ax

    def plot_site_scatter(self, df_results, site_name, save_plots=True):
        """
        Create scatter plots for each AERONET site showing TROPOMI vs AERONET AOD
        
        Parameters:
        -----------
        df_results : pandas.DataFrame
            DataFrame containing collocation results
        site_name : str
            Name of the AERONET site to plot
        target_wavelength : int
            Target wavelength for AERONET AOD (e.g., 380, 500, 675 nm)
        save_plots : bool
            Whether to save plots to file
        """
        if df_results is None or len(df_results) == 0:
            print("No data available for scatter plots")
            return
        
        
        try:
            # Filter data for the specific site
            site_data = df_results[df_results['site_name'] == site_name].copy()
            
            if len(site_data) == 0:
                print(f"  - No data for {site_name}")
                return None, None  # Return consistent tuple format 
            
            # Convert datetime if needed
            if isinstance(site_data['aeronet_datetime'].iloc[0], str):
                site_data['aeronet_datetime'] = pd.to_datetime(site_data['aeronet_datetime'])
            
            # Get site information
            site_info = self._get_site_info_string(site_data, site_name)
            
            # Get time range
            start_date = site_data['aeronet_datetime'].min()
            end_date = site_data['aeronet_datetime'].max()
            time_range = f"{start_date.strftime('%Y-%m-%dT%H')} - {end_date.strftime('%Y-%m-%dT%H')}"
            
            # Get wavelength information
            tropomi_wavelength = site_data['tropomi_wavelength'].iloc[0]
            
            # Create the scatter plot
            fig, ax = plt.subplots(figsize=(10, 10))
            
            # Main scatter plot
            scatter = ax.scatter(site_data['aeronet_aod_mean'], site_data['tropomi_aod_mean'],
            #                   c=range(len(site_data)), cmap='viridis', 
            #                   s=60, alpha=0.7, edgecolors='black', linewidth=0.5)
                                s=60, alpha=0.7, edgecolors='black', linewidth=0.5,
                                color='gray')  # Simple single color instead of temporal mapping
            
            # Add error bars for TROPOMI if available
            if 'tropomi_aod_std' in site_data.columns:
                ax.errorbar(site_data['aeronet_aod_mean'], site_data['tropomi_aod_mean'],
                           yerr=site_data['tropomi_aod_std'], fmt='none', 
                           color='gray', alpha=0.5, capsize=3, elinewidth=0.8)
 
            # Add error bars for AERONET if available
            if 'aeronet_aod_std' in site_data.columns:
                ax.errorbar(site_data['aeronet_aod_mean'], site_data['tropomi_aod_mean'],
                           xerr=site_data['aeronet_aod_std'], fmt='none', 
                           color='gray', alpha=0.5, capsize=3, elinewidth=0.8)
                
            # Calculate limits for 1:1 line
            min_val = min(site_data['aeronet_aod_mean'].min(), site_data['tropomi_aod_mean'].min())
            max_val = max(site_data['aeronet_aod_mean'].max(), site_data['tropomi_aod_mean'].max())
            
            # Add some padding
            padding = (max_val - min_val) * 0.1
            plot_min = max(0, min_val - padding)  # Don't go below 0 for AOD
            plot_max = max_val + padding
            
            # 1:1 line
            ax.plot([plot_min, plot_max], [plot_min, plot_max], 
                   'r--', linewidth=2, alpha=0.8, label='1:1 Line')
            
            # Calculate and plot regression line
            from scipy import stats
            slope, intercept, r_value, p_value, std_err = stats.linregress(
                site_data['aeronet_aod_mean'], site_data['tropomi_aod_mean'])
            
            regression_x = np.array([plot_min, plot_max])
            regression_y = slope * regression_x + intercept
            ax.plot(regression_x, regression_y, 'b-', linewidth=2, alpha=0.8,
                   label=f'Regression (y = {slope:.3f}x + {intercept:.3f})')
            
            # Set equal aspect ratio and limits
            ax.set_xlim(plot_min, plot_max)
            ax.set_ylim(plot_min, plot_max)
            ax.set_aspect('equal')
            
            # Labels and formatting
            ax.set_xlabel(f'AERONET AOD ({self.target_wavelength}nm)', fontsize=12, fontweight='bold')
            ax.set_ylabel(f'TROPOMI AOD ({tropomi_wavelength}nm)', fontsize=12, fontweight='bold')
            
            # Get matching criteria information
            matching_info = f"{self.spatial_radius_km}km and ±{self.temporal_window_hours}hr"

            # Title with site info and time range
            title = f'AERONET AOD Scatter Plot - Matching TROPOMI within {matching_info} \nSite:{site_name}{site_info}\nTime_range:{time_range}'
            ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
            
            # Grid and legend
            ax.grid(True, alpha=0.3)
            ax.legend(loc='upper left', fontsize=10)
            
            # Add colorbar for temporal information - You can see if the TROPOMI-AERONET agreement changes over time
            #cbar = plt.colorbar(scatter, ax=ax, shrink=0.8)
            #cbar.set_label('Temporal Order', fontsize=10)
            
            # Calculate and display statistics
            correlation = np.corrcoef(site_data['aeronet_aod_mean'], site_data['tropomi_aod_mean'])[0,1]
            rmse = np.sqrt(np.mean((site_data['aeronet_aod_mean'] - site_data['tropomi_aod_mean'])**2))
            bias = np.mean(site_data['tropomi_aod_mean'] - site_data['aeronet_aod_mean'])
            mae = np.mean(np.abs(site_data['aeronet_aod_mean'] - site_data['tropomi_aod_mean']))
            n_points = len(site_data)
            
            # Statistics text box
            stats_text = (f'N = {n_points}\n'
                         f'R = {correlation:.3f}\n'
                         f'RMSE = {rmse:.3f}\n'
                         f'Bias = {bias:.3f}\n'
                         f'MAE = {mae:.3f}\n'
                         f'Slope = {slope:.3f}')
            
            ax.text(0.98, 0.02, stats_text, transform=ax.transAxes, fontsize=10,
                   verticalalignment='bottom', horizontalalignment='right',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
            
            plt.tight_layout()
            
            # Save plot if requested
            if save_plots:
                plot_filename = f"scatter_{site_name.replace(' ', '_').replace('-', '_')}.png"
                plot_path = os.path.join(self.output_dir, 'plots', plot_filename)
                plt.savefig(plot_path, dpi=300, bbox_inches='tight')
                print(f"     Scatter plot saved for {site_name}: {plot_path}")
            
            plt.show()
            
        except Exception as e:
            print(f"  - Error creating scatter plot for {site_name}: {str(e)}")
    
    def _get_site_info_string(self, site_data, site_name):
        """
        Helper function to format site information string
        
        Parameters:
        -----------
        site_data : pandas.DataFrame
            Site-specific data
        site_name : str
            Name of the site
            
        Returns:
        --------
        str
            Formatted site information string
        """
        site_info = ""
        try:
            # Try to get from aeronet_sites first
            if hasattr(self, 'aeronet_sites') and site_name in self.aeronet_sites:
                site_coords = self.aeronet_sites[site_name]
                
                if isinstance(site_coords, dict):
                    lat = site_coords.get('lat')
                    lon = site_coords.get('lon')
                    elevation = site_coords.get('elevation')
                    
                    if lat is not None and lon is not None:
                        site_info = f" (lat: {lat:.2f}°, lon: {lon:.2f}°"
                        if elevation is not None:
                            site_info += f", Alt: {elevation:.1f}m)"
                        else:
                            site_info += ")"
            
            # Fallback to DataFrame coordinates if aeronet_sites didn't work
            if not site_info and len(site_data) > 0:
                if 'site_lat' in site_data.columns and 'site_lon' in site_data.columns:
                    lat = site_data['site_lat'].iloc[0]
                    lon = site_data['site_lon'].iloc[0]
                    site_info = f" (lat: {lat:.2f}°, lon: {lon:.2f}°)"
                    
        except Exception as e:
            print(f"Warning: Could not format site info for {site_name}: {str(e)}")
        
        return site_info

    def plot_validation_scatter(self, df_results, save_plots=True):
        """
        Create a single scatter plot showing TROPOMI vs AERONET AOD for all collocated sites
        
        Parameters:
        -----------
        df_results : pandas.DataFrame
            DataFrame containing collocation results from all sites
        target_wavelength : int
            Target wavelength for AERONET AOD (e.g., 380, 500, 675 nm)
        save_plots : bool
            Whether to save plots to file
        """
        if df_results is None or len(df_results) == 0:
            print("No data available for validation scatter plot")
            return
        
        print(f"\nGenerating validation scatter plot for all {len(df_results['site_name'].unique())} sites...")
        
        try:
            # Convert datetime if needed
            if isinstance(df_results['aeronet_datetime'].iloc[0], str):
                df_results['aeronet_datetime'] = pd.to_datetime(df_results['aeronet_datetime'])
            
            # Get time range for all data
            start_date = df_results['aeronet_datetime'].min()
            end_date = df_results['aeronet_datetime'].max()
            time_range = f"{start_date.strftime('%Y-%m-%d')} - {end_date.strftime('%Y-%m-%d')}"
            
            # Get wavelength information
            tropomi_wavelength = df_results['tropomi_wavelength'].iloc[0]
            
            # Get unique sites for color coding
            unique_sites = df_results['site_name'].unique()
            #colors = plt.cm.tab10(np.linspace(0, 1, len(unique_sites)))
            colors = plt.cm.Set1(np.linspace(0, 1, len(unique_sites))) # 9 very distinct colors
            site_color_map = dict(zip(unique_sites, colors))
            
            # Create the scatter plot
            fig, ax = plt.subplots(figsize=(12, 10))
            
            # Plot each site with different colors
            for i, site_name in enumerate(unique_sites):
                site_data = df_results[df_results['site_name'] == site_name]
                
                scatter = ax.scatter(site_data['aeronet_aod_mean'], site_data['tropomi_aod_mean'],
                                   s=60, alpha=0.7, edgecolors='black', linewidth=0.5,
                                   color=colors[i], label=f'{site_name} (n={len(site_data)})')
                
                # # Add error bars for TROPOMI if available
                # if 'tropomi_aod_std' in site_data.columns:
                #     ax.errorbar(site_data['aeronet_aod_mean'], site_data['tropomi_aod_mean'],
                #                yerr=site_data['tropomi_aod_std'], fmt='none', 
                #                color=colors[i], alpha=0.3, capsize=2, elinewidth=0.5)
                # # Add error bars for AERONET if available
                # if 'aeronet_aod_std' in site_data.columns:
                #     ax.errorbar(site_data['aeronet_aod_mean'], site_data['tropomi_aod_mean'],
                #                xerr=site_data['aeronet_aod_std'], fmt='none', 
                #                color=colors[i], alpha=0.3, capsize=2, elinewidth=0.5)
                    
            # Calculate limits for 1:1 line
            min_val = min(df_results['aeronet_aod_mean'].min(), df_results['tropomi_aod_mean'].min())
            max_val = max(df_results['aeronet_aod_mean'].max(), df_results['tropomi_aod_mean'].max())
            
            # Add 1:1 line
            ax.plot([min_val, max_val], [min_val, max_val], 'r--', 
                   linewidth=2, alpha=0.8, label='1:1 Line')
            
            # Calculate and display statistics for all data
            correlation = df_results['aeronet_aod_mean'].corr(df_results['tropomi_aod_mean'])
            rmse = np.sqrt(np.mean((df_results['aeronet_aod_mean'] - df_results['tropomi_aod_mean'])**2))
            bias = np.mean(df_results['tropomi_aod_mean'] - df_results['aeronet_aod_mean'])
            n_points = len(df_results)
            
            # Add statistics text
            stats_text = (f'All Sites Statistics:\n'
                         f'N = {n_points}\n'
                         f'R = {correlation:.3f}\n'
                         f'RMSE = {rmse:.3f}\n'
                         f'Bias = {bias:.3f}')
            
            ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
                   verticalalignment='top', bbox=dict(boxstyle='round', 
                   facecolor='white', alpha=0.8), fontsize=12)
            
            # Labels and title
            ax.set_xlabel(f'AERONET AOD @ {self.target_wavelength} nm', fontsize=14, fontweight='bold')
            ax.set_ylabel(f'TROPOMI AOD @ {tropomi_wavelength} nm', fontsize=14, fontweight='bold')
            ax.set_title(f'TROPOMI vs AERONET AOD Validation\nAll Sites Combined\n{time_range}', 
                        fontsize=16, fontweight='bold')
            
            # Set equal aspect ratio and limits
            ax.set_aspect('equal', adjustable='box')
            buffer = (max_val - min_val) * 0.05
            ax.set_xlim(min_val - buffer, max_val + buffer)
            ax.set_ylim(min_val - buffer, max_val + buffer)
            
            # Add grid
            ax.grid(True, alpha=0.3)
            
            # Add legend
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
            
            plt.tight_layout()
            
            if save_plots:
                if not os.path.exists(f"{self.output_dir}/plots"):
                    os.makedirs(f"{self.output_dir}/plots")
                
                filename = f"{self.output_dir}/plots/validation_scatter_all_sites_{self.target_wavelength}nm.png"
                plt.savefig(filename, dpi=300, bbox_inches='tight')
                print(f"  - Saved validation scatter plot: {filename}")
            
            plt.show()
            
        except Exception as e:
            print(f"  - Error creating validation scatter plot: {str(e)}")
        
