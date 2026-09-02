# spec/spec_helper.rb
require 'rack/test'
require 'rspec'
require_relative '../lib/server.rb'  # <-- load the proper app entry point

ENV['RACK_ENV'] = 'test'

module SinatraHelper
  def app
    Server.new
  end
end

RSpec.configure do |config|
  config.include SinatraHelper
  config.include Rack::Test::Methods
end
