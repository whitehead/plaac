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

### Locally

1. Create a ```logs``` directory:```mkdir logs```
2. Use the included **dev_server** script to launch a  development webserver at [http://localhost:4567](http://localhost:4567)
   (note that **shotgun_server** is now obsolete)

```bash
bin/dev_server
```

### In Apache/Nginx

The server is a rack compatable application. Please consult the official guides for your webserver, e..g Ngnix or Apache.

### On a server

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

```bash
sudo apt update && sudo apt install -y git
```

#### Initial setup

Clone the PLAAC repository on the server:

```bash
git clone https://github.com/whitehead/plaac.git
cd plaac/web
```

Run the server setup script:

```bash
./deploy/setup-plaac-server.sh
```

The script, installs podman builds the PLAAC web container from the current checkout, installs the Podman/Quadlet service, enables the service to run without an interactive login, and starts the application.

Check the service:

```bash
systemctl --user status plaac.service
```

Check the container:

```bash
podman ps
```

### Accessing the application

By default, the application listens on TCP port `4567`.

From the server itself:

```bash
curl http://127.0.0.1:4567/
```

A successful response should return HTTP 200. To find the server's IP address on Debian:

```bash
hostname -I
```

For example:

```text
192.168.122.42
```

From another computer on the same network, access the application at:

```text
http://192.168.122.42:4567/
```

For a VPS, `hostname -I` may show the server's private address rather than its public Internet address. The VPS provider's control panel should be used to determine the public IP address when necessary.

If the server is protected by a firewall, TCP port 4567 must be allowed. On a Debian server using UFW:

```bash
sudo ufw allow 4567/tcp
```

Cloud VPS providers may also have a separate network firewall or security-group configuration. If so, TCP port 4567 must be allowed there as well.

Port 4567 is suitable for initial testing. For a public production deployment, the application should normally be placed behind a reverse proxy and HTTPS should be used instead of exposing the application directly.

### Updating an existing deployment

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
