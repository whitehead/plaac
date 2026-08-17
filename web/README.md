# Developer installation for PLAAC web application

## Requirements

#### Java
See [java.com](http://www.java.com/en/) or reference your OS's installation instructions.

#### Rscript & dependencies

On Debian:
```bash
sudo apt-get install littler r-cran-cairodevice
```
On Fedora:
```bash
sudo yum install R-core
```

#### Ruby >= 3.3.10 ([with rvm](https://rvm.io/))
```bash
rvm install 3.3.10
rvm use 3.3.10
```

#### Ruby Dependencies ([bundler](http://bundler.io/))
```bash
gem install bundler
bundle install
```

## Deploying the web application

### Local deployment

1. Create a ```logs``` directory:```mkdir logs```
2. Ensure the `plaac.jar` is built, following the instructions in [`cli/README.md`](../cli/README.md)
3. Copy the `plaac.jar` to the right place:
```
cp ../cli/target/plaac.jar bin/plaac.jar
```
4. Use the included **dev_server** script to launch a  development webserver at [http://localhost:4567](http://localhost:4567)
   (note that **shotgun_server** is now obsolete)

```bash
bin/dev_server
```

### In Apache/Nginx

The server is a rack compatable application. Please consult the official guides for your webserver, e..g Ngnix or Apache.

### Container-based deployment

The web application can be deployed as a Podman container on a Linux
server. The container includes the Ruby application and its
dependencies, the Java runtime required to run `plaac.jar`, and R for
generating plots. The server therefore does not need Ruby, Java, or R
installed directly.

The deployment uses Podman and systemd/Quadlet. The container image is built locally from the PLAAC source checkout, so access to a container registry is not required.

#### Server requirements

The server requires:

* Linux with systemd
* Git
* A user account with permission to run Podman
* Internet access to clone the PLAAC repository and build the container image

The following instructions assume a Debian distribution:


#### Initial setup

Install Git:

```bash
sudo apt update && sudo apt install -y git
```

Clone the PLAAC repository on the server:

```bash
git clone https://github.com/whitehead/plaac.git
cd plaac/web
```

Run the server setup script:

```bash
./deploy/setup-plaac-server.sh
```

The script installs `podman`, builds the PLAAC web container from the current checkout, installs the Podman/Quadlet service, enables the service to run without an interactive login, and starts the application.

Check the service:

```bash
systemctl --user status plaac.service
```

Check the container:

```bash
podman ps
```

#### Accessing the application

By default, the application listens on TCP port `4567`.

From the server itself:

```bash
curl -I http://127.0.0.1:4567/
```

A successful response should return HTTP 200. 

The `setup-plaac-server.sh` also configures Caddy as a reverse proxy, so the application can be accessed externally through HTTP without exposing its internal port directly.

By default, access the server using its public IP address:

```text
http://<server-ip>/
```

The server's IP address can be found with:

```bash
hostname -I
```

For a public deployment, a domain name and access via HTTPS can be configured by setting the `DOMAIN` environment variable:

```bash
DOMAIN=plaac.example.com ./deploy/setup-plaac-server.sh
```

Caddy will then configure HTTPS automatically:

```text
https://plaac.example.com/
```

The domain's DNS record must point to the server's public IP address. If the server's IP address changes, only the DNS record needs to be updated.

If `DOMAIN` is not set, the server remains available through its IP address over HTTP.

> The application listens internally on TCP port `4567`. Caddy proxies requests to this port and listens publicly on port `80`. When a domain is configured, Caddy also listens on port `443` and manages the TLS certificate.  Cloud VPS providers may also have a separate network firewall or security-group configuration. If so, TCP port 4567 must be allowed there as well.

The Caddy HTTP configuration does not contain the server's IP address, so it continues to work if the public IP changes.

#### Updating an existing deployment

To update the application, update the source checkout and run the setup script again:

```bash
cd ~/plaac
git pull
cd web
./deploy/setup-plaac-server.sh
```

The script rebuilds the container from the updated source and restarts the service.

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
