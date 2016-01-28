# Developer installation for PLAAC web application

## Requirements


#### Java
See [java.com](http://www.java.com/en/) or reference your OS's installation instructions.

#### Rscript & dependencies
On Debian:
```
    $ sudo apt-get install littler r-cran-cairodevice
```
On Fedora:
```
    $ sudo yum install R-core
```

#### Ruby >= 1.9.3 ([with rvm](https://rvm.io/))
```
    $ rvm install 1.9.3
    $ rvm use 1.9.3
```

#### Ruby Dependencies ([bundler](http://bundler.io/))
```
    $ gem install bundler
    $ bundle install
```
-------------------------------------------------

## Running The Server


#### Locally

1. create a ```logs``` directory:```mkdir logs```
2. Use the included **shotgun_server** script to launch a  development webserver at [http://localhost:4567](http://localhost:4567)

```
    $ bin/shotgun_server
```

#### In Apache/Nginx

The server is a rack compatable application. Please consult the official guides for your webserver.

* [Phusion Passenger Nginx Guide](http://www.modrails.com/documentation/Users%20guide%20Apache.html)
* [Phusion Passenger Apache Guide](http://www.modrails.com/documentation/Users%20guide%20Apache.html)

-------------------------------------------------

## Contributing/Hacking

Pull requests are welcome.

#### Code Layout

This is a [sinatra](http://www.sinatrarb.com/) app. Start by looking at  ```lib/server.rb``` 

* ```lib/server.rb``` - Where all the url controller methods live, along with the building of the commandline options for the plaac.jar

* ```views/*.haml``` - All the webpage templates. See [haml docs](http://haml.info/) for more info.

* ```views/js/*.coffee``` - All the javascript behavior, written in [coffeescript](http://coffeescript.org/).

#### Application Flow

Designed to handle the submission of FASTA files
to a sequence analysis application called plaac


**flow**
```
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
