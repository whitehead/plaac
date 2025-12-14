require 'sinatra/extension'
require 'coffee-script'
require 'haml'
require 'sass'

module Assets
  extend Sinatra::Extension

  Tilt.register Tilt::CoffeeScriptTemplate, 'coffee'
  
  @@base = File.expand_path('../..',__FILE__)

  set :public_folder, @@base+"/public"
  set :views,         @@base+"/views"

  get '/sass/:script.css' do
    clean! params[:script]
    sass :"css/#{params[:script]}"
  end

  get '/coffee/:script.js' do
    clean! params[:script]
    #coffee :"js/#{params[:script]}"

    file_path = File.join(settings.views, 'js', "#{params[:script]}.coffee")
    halt 404, "Not found" unless File.exist?(file_path)
    content_type 'application/javascript'
    CoffeeScript.compile(File.read(file_path))
  end

  get '/js/:script' do
    clean! params[:script]
    catfile "views/js/#{params[:script]}"
  end

  get '/:haml.html' do
    clean! params[:haml]
    haml :"views/#{params[:haml]}"
  end

  helpers do
    def catfile(filename)
      return nil if !File.file?(filename)
      File.open(filename,'r') {|f| f.read}
    end

    def clean!(filename)
      return nil if filename.nil? || filename.empty?
      filename.gsub!(/(\.\.|[\/\?\*])/, '')
    end
  end
end
