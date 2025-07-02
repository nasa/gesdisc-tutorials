# Custom GES DISC Anaconda Environments
#### Description
In order to mitigate errors, GES DISC provides customized Anaconda YAML environment files, which users can install before running select Python tutorials. These YAMLs are organized by services, such as OPeNDAP, and can be installed directly or by downloading from the gesdisc-tutorials/environments folder. The instructions for installing each individual environment are described here, and will be updated as more environment files are added or changed.

#### 1. Installing Anaconda
An Anaconda installation is required for installing custom environment files. We recommend installing `miniconda` for lightweight, terminal-based environment management, or the standard Anaconda Distribution, which includes a GUI (Anaconda Navigator) for managing environment files. `miniconda` installation steps can be accessed [here](https://www.anaconda.com/docs/getting-started/miniconda/install#quickstart-install-instructions), and Anaconda Distribution steps can be accessed [here](https://www.anaconda.com/docs/getting-started/anaconda/install), available for Windows/macOS/Linux.

#### 2. Installing `nasa-gesdisc.yml`
The `nasa-gesdisc.yml` environment contains a Python 3.10 installation that works with all of our how-tos *except* for our [Cloud OPeNDAP guides](https://disc.gsfc.nasa.gov/information/tools?keywords=cloud%20opendap&title=OPeNDAP%20In%20The%20Cloud). If you wish to access Cloud OPeNDAP, we recommend following the steps for installing `opendap.yml` outlined in step 3.

- To install via a command line: 
  - Open a terminal, ensuring that the `conda` can be called.
  - Enter: `conda env create -f https://raw.githubusercontent.com/nasa/gesdisc-tutorials/refs/heads/main/environments/nasa-gesdisc.yml`.
  - Activate the environment by entering: `conda activate nasa-gesdisc`.

- To install via Anaconda Navigator: 
  - Download `nasa-gesdisc.yml` to your computer.
  - Open Anaconda Navigator, click `Environments` -> `Import`
  - Select `Import from: Local drive` and select the downloaded `nasa-gesdisc.yml` file. `New Environment Name` will automatically fill using the name from the file, but you can change this if desired.
  - Click the triangular "play" button to the right of the `nasa-gesdisc` environment tab, and click the "Open in Terminal" button to activate the environment in a terminal session

#### 3. Installing `opendap.yml`
The `opendap.yml` environment contains a Python 3.10 installation that *only* works with select OPeNDAP how-tos, which includes the `Pydap` library. If you are not accessing an OPeNDAP service with Python, we recommend installing the `nasa-gesdisc` environment, outlined in step 2.

- To install via a command line: 
  - Open a terminal, ensuring that the `conda` can be called.
  - Enter: `conda env create -f https://raw.githubusercontent.com/nasa/gesdisc-tutorials/refs/heads/main/environments/opendap.yml`.
  - Activate the environment by entering: `conda activate nasa-gesdisc-opendap`.

- To install via Anaconda Navigator: 
  - Download `opendap.yml` to your computer.
  - Open Anaconda Navigator, click `Environments` -> `Import`
  - Select `Import from: Local drive` and select the downloaded `opendap.yml` file. `New Environment Name` will automatically fill using the name from the file, but you can change this if desired.
  - Click the triangular "play" button to the right of the `opendap` environment tab, and click the "Open in Terminal" button to activate the environment in a terminal session.