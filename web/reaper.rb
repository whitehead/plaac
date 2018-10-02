#!/usr/bin/env ruby
require 'date'
require 'fileutils'
require 'logger'
log = Logger.new('reaper.log')
if Dir.exists? "/tmp/clusty"
  Dir["/tmp/clusty/*"].each do |dir|
    next unless File.directory? dir

    stats = File::Stat.new(dir)
    age_threshold_seconds = 60*60*24
    age_threshold = Time.now - age_threshold_seconds
    if stats.ctime < age_threshold
      log.info "deleting #{dir} ctime=#{stats.ctime}"
      FileUtils.rm_rf(dir)
    end
  end
end
