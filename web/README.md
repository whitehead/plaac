# PLAAC web application

The PLAAC web application can be run in several ways, depending on what you need to do:

- **Run a prebuilt container locally** - the simplest option for running and testing the web application using docker. No application dependencies need to be installed on the host.
- **Deploy a container on a server** - intended primarily for website administrators who need to deploy PLAAC on an external server or hosting environment.
- **Install and run the application manually**  intended for development or environments where container deployment is not appropriate. This requires installing the application's dependencies directly.

The recommended approach for most users who just want to run the web application locally is the prebuilt container.

## Running the prebuilt container locally

A prebuilt PLAAC web container can be run locally without installing
Ruby, Java, R, or the other application dependencies required to run
`plaac.jar`.  This is useful for development and testing on macOS,
Windows, or Linux systems with [Docker](https://www.docker.com/) or
[Podman](https://podman.io/) installed.

The prebuilt containers are published for both `linux/amd64` and
`linux/arm64` and therefore run natively on Apple Silicon Macs.

For example, using Docker:

```bash
docker run --rm -p 4567:4567 ghcr.io/<username>/plaac-web:latest
```

Or using Podman:

```bash
podman run --rm -p 4567:4567 ghcr.io/<username>/plaac-web:latest
```

Replace `<username>` with the GitHub username or organization that publishes the PLAAC container image.

On macOS, the recommended installation path for Podman is via the [dedicated installer](https://podman.io/docs/installation), they do not recommend using [Homebrew](https://brew.sh/).  [podman desktop](https://podman-desktop.io/downloads), a GUI version, is also available.

```bash
podman machine init
podman machine start
```

Then run the container as above.

Access the application at [http://localhost:4567](http://localhost:4567).

The `-p 4567:4567` option publishes the application's port from the container to the local machine. No reverse proxy is required for local use.

On macOS and Windows, Docker or Podman provides the Linux environment required to run the container; no separate Debian installation or manually managed server VM is required.

## Deploy a container on a server

The web application can be deployed as a persistent service on a
Linux server using Podman. The deployment uses systemd/Quadlet and
can optionally configure Caddy as a reverse proxy for public access
and HTTPS.

The container image can either be pulled from a container registry or
built locally from a PLAAC source checkout, if no container image is
supplied. The server itself does not need Ruby, Java, or R installed
directly.

### Server requirements

The server requires:

* Linux with systemd
* Git
* A user account with permission to run Podman
* Internet access to clone the PLAAC repository and build the container image

The following instructions assume a Debian distribution:

### Initial setup

#### 1. Install Git:

```bash
sudo apt update && sudo apt install -y git
```

#### 2. Clone the PLAAC repository on the server

```bash
git clone https://github.com/whitehead/plaac.git
cd plaac/web
```

#### 3. Run the server setup script:

```bash
./deploy/setup-plaac-server.sh
```

The script installs `podman`, builds or uses the specified PLAAC web container, installs the Podman/Quadlet service, enables the service to run without an interactive login, and starts the application.

Optional arguments/environment variables:

* **Supply container image**.  By default the script builds the PLAAC web container from the current checkout. Alternatively, a prebuilt container image can be provided as a command-line argument:

  ```bash
  ./deploy/setup-plaac-server.sh <container>
  ```

  For example, `<container>` can be an image reference supported by Podman, in the format of `ghcr.io/<username>/plaac-web:latest`.

* **Supply domain name**. You can also optionally supply a domain name to allow access via HTTPS, suitable for public deployment. This is configured by setting the `DOMAIN` environment variable:

  ```bash
  DOMAIN=plaac.example.com ./deploy/setup-plaac-server.sh
  ```

#### 4. Check the service:

```bash
systemctl --user status plaac.service
```

#### 5. Check the container:

```bash
podman ps
```

### Accessing the application

By default, the application listens on TCP port `4567`.

To check initial setup, from the server itself run:

```bash
curl -I http://127.0.0.1:4567/
```

A successful response should return `HTTP 200`. 

The `setup-plaac-server.sh` also configures Caddy as a reverse proxy, so the application can be accessed externally through HTTP without exposing its internal port directly.

By default, access the server using its public IP address:

```text
http://<server-ip>/
```

The server's IP address can be found with:

```bash
hostname -I
```

If you had set the `DOMAIN` environment variable during the setup above, Caddy will then configure HTTPS automatically with whatever you had used, e.g.:

```text
https://plaac.example.com/
```

The domain's DNS record must point to the server's public IP address. If the server's IP address changes, only the DNS record needs to be updated.  If `DOMAIN` is not set, the server remains available through its IP address over HTTP.

> The application listens internally on TCP port `4567`. Caddy proxies requests to this port and listens publicly on port `80`. When a domain is configured, Caddy also listens on port `443` and manages the TLS certificate.  Cloud VPS providers may also have a separate network firewall or security-group configuration. If so, TCP port 4567 must be allowed there as well.

The Caddy HTTP configuration does not contain the server's IP address, so it continues to work if the public IP changes.

### Updating an existing deployment

To update the application, update the source checkout and run the setup script again:

```bash
cd ~/plaac
git pull
cd web
./deploy/setup-plaac-server.sh
```

The script rebuilds the container from the updated source and restarts the service. If you had supplied any optional container or domain name options in the original run, you should supply them again, otherwise the script will fallback to defaults. 

## Install and run the application manually

### Requirements

The PLAAC web application requires the following software and dependencies:

- **Java**, required to run the PLAAC Java application.
- **Rscript**, with the R packages required by PLAAC, including Cairo support.
- **Ruby >= 3.3.10**.
- **Bundler** and the Ruby dependencies specified by the application.
- A built **`plaac.jar`** from the PLAAC CLI component.

The exact installation method for these dependencies depends on the deployment environment. The local installation instructions below provide an example for a development environment.

## Local deployment of web application

The following steps set up the PLAAC web application for running on a local machine, directly from a git checkout.

### 1. Install the required system dependencies

Install R and graphviz (needed for `dot` which generates HMM probability image), on Debian-based systems:

```bash
sudo apt-get install littler r-cran-cairodevice graphviz
```

On Fedora:

```bash
sudo yum install R-core graphviz
```

Install Java using [java.com](http://www.java.com/en/) or your
operating system's package manager and installation instructions. You
will need both the JRE (for running the JAR) and the JDK (to compile
it in the `build_plaac.sh` step).

On Debian:

```bash
sudo apt-get install openjdk-17-jre-headless openjdk-17-jdk-headless
```

On Fedora:

```bash
sudo yum install java-17-openjdk-headless java-17-openjdk-devel
```

### 2. Install Ruby

PLAAC requires Ruby >= 3.3.10, typically via your distribution packages:

On Debian:

```bash
sudo apt-get install ruby
```

On Fedora:

```bash
sudo yum install ruby rubygems
```

If using [RVM](https://rvm.io/):

```bash
rvm install 3.3.10
rvm use 3.3.10
```

### 3. Install the Ruby dependencies

Install Bundler and the application's Ruby dependencies:

```bash
cd web
gem install bundler
bundle install
cd ..
```

### 4. Create the logs directory

```bash
mkdir logs
```

### 5. Build the PLAAC JAR

Ensure that `plaac.jar` is built (detailed instructions are in [`cli/README.md`](../cli/README.md)) and copy the resulting JAR and other generated files into the web application directory:

```bash
cd cli/
./build_plaac.sh
cp target/plaac.jar ../web/bin/
cp target/_plaac_headers.haml ../web/views/
cp target/hmm_default.png ../web/public/
cd ..
```
### 6. Start the development server

Use the included **`dev_server`** script to launch the development webserver at [http://localhost:4567](http://localhost:4567):

```bash
cd web/
./bin/dev_server
```

The `shotgun_server` script is obsolete and should not be used.

> **Apache/Nginx** The server is a rack compatable application. Please consult the official guides for your webserver, e..g Ngnix or Apache.

## Contributing/Hacking

Pull requests are welcome.

#### Code Layout

This is a [sinatra](http://www.sinatrarb.com/) app. Start by looking at  ```lib/server.rb``` 

* ```lib/server.rb``` - Where all the url controller methods live, along with the building of the commandline options for the plaac.jar

* ```views/*.haml``` - All the webpage templates. See [haml docs](http://haml.info/) for more info.

* ```views/js/*.coffee``` - All the javascript behavior, written in [coffeescript](http://coffeescript.org/).

#### Application flow

Designed to handle the submission of FASTA files to the command-line `plaac.jar` program.

**flow**
```text
get '/' - input form and options
post '/find' - runs plaac.jar
get '/candidates/:token' - lists candidates, visualization input form
post '/visualize' - runs plaac.jar for details and rscripts to visualize
get '/visualize/:token' - visualization output

# Downloadable content
get '/candidates/:token/tsv' 
get '/visualize/:token/pdf'
get '/visualize/:token/images/:id.png' do
get '/visualize/:token/tsv' do

get '/details' - algorithm description and links.
```
