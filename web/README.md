Developer installation for PLAAC web application
================================================

Requirements
------------

* Ruby 1.9.3 https://rvm.io/
    > rvm install 1.9.3
    > rvm use 1.9.3

* Dependencies (bundler)
    > gem install bundler
    > bundle install

Infrastructure
--------------

* bin/sinatra_server
 *  a webrick server (same as ruby lib/server.rb)

* bin/shotgun_server
 *  a development server that reloads all the 
    code every request (slow but useful)

* config.ru
 *  makes the application rackable 
    (e.g. you can drop it into passenger)

* lib/server.rb
 *  Start here.

Application
-----------

Designed to handle the submission of FASTA files
to a sequence analysis application called plaac


* flow:
 *  / - input form
 *  /run - handles post from the input form
         runs plaac with given options
 *  /sequence(.*).fasta - the format of the results files
