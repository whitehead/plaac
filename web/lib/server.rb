# Encoding: utf-8
require 'bundler/setup'
require 'sinatra/base'
require 'sinatra/cookies'
require 'fileutils'
require 'logger'
require 'bio'
require 'yaml'
require 'uuid'
require 'better_errors'
require 'benchmark'

$:.push File.expand_path("..",__FILE__)
require 'assets'
require 'job'

#if RUBY_VERSION =~ /1.9/ # assuming you're running Ruby ~1.9
#  Encoding.default_external = Encoding::UTF_8
#  Encoding.default_internal = Encoding::UTF_8
#end

# ensure required directories
%w[tmp logs].each do |dir|
  dir = File.expand_path("../../"+dir,__FILE__)
  FileUtils.mkdir(dir) unless File.directory?(dir)
end

ENVIRONMENT = ENV["RAILS_ENV"] || "development"
$version = `git log | head -n 1`.split(/\s/)[1][0,10] rescue "?"

class Server < Sinatra::Base
  helpers Sinatra::Cookies
  register Assets

  enable :sessions
  set :session_secret, 'marshmallow'
  set :bind, '0.0.0.0'

  configure :development do
    use BetterErrors::Middleware
    BetterErrors.application_root = File.dirname(__FILE__)
  end

  $config = YAML.load_file(File.expand_path("../../config.yaml",__FILE__))

  @@base = File.expand_path('../..',__FILE__)

  @@log = Logger.new(@@base+"/logs/plaac.log")

  @@input_fasta="input_sequence.fasta"
  @@prion_candidates="prion_candidates.tsv"
  @@prion_candidates_details="prion_candidates_details.tsv"
  @@prion_candidates_details_pdf="prion_candidates.pdf"

  @@background_meta = "bg_freqs/background_metadata.txt"

  bg_freq_filename = ->(id){ id.nil? ? nil : "bg_freqs_#{id}.txt" }
  bg_freq_local_filename = ->(id){ id.nil? ? nil : "bg_freqs/#{bg_freq_filename[id]}" }

  before do
    # Strip the last / from the path
    request.env['PATH_INFO'].gsub!(/\/$/, '')
  end

  #--------------------------------------------------------------------
  # Pages

  get '/' do
    @bg_freqs = read_bg_freqs()
    haml :index
  end

  #--------------------------------------------------------------------
  # Actions

  # upload a fasta file for processing
  post '/find' do
    @@log.info("POST /run")
    @bg_freqs = read_bg_freqs()
    # validation
    unless params[:len] && params[:len].to_i > 0
      @error = "Core length must be a positive integer"
      @@log.info(@error)
      return haml(:index)
    end

    session[:core_len] = core_len = params[:len].to_i
    session[:alpha] = alpha = params[:alpha].to_f / 100.0 rescue 0.0

    job = Job.new
    working_directory = job.working_directory
    output_filename = File.join(working_directory,@@prion_candidates)
    filename = File.join(working_directory,@@input_fasta)

    job.add_file(path("bin/plaac_plot.r"), working_directory) # r library
    job.add_file(path("bin/plaac_plot_util.r"), working_directory) # r library
    job.add_file(path("bin/plaac.jar"), working_directory) # binary

    # Create input file on disk.
    if params[:sequence] && !params[:sequence].empty?
      # Accept String (pasted)
      write_string(filename, add_gene_id(params[:sequence]))
    elsif params[:file] &&
        (tmpfile = params[:file][:tempfile]) &&
        (name = params[:file][:filename])
      # Accept File (uploaded)
      write_file(filename, tmpfile)
    else
      @error = "No sequence."
      @@log.info(@error)
      return haml(:index)
    end

    @@log.debug "validating fasta #{filename}"
    valid_fasta, msg = valid_fasta_file(filename)
    @@log.debug "valid?: #{valid_fasta}, #{msg}"
    if !valid_fasta
      @error = msg
      return haml(:index)
    end

    gene_count = `egrep '^>' #{filename} |wc -l`.to_i

    @@log.info "gene_count: #{gene_count}"

    bg_freqs = read_bg_freqs() # ignore warning, bg_freqs used by filename lambdas
    bgfreq_filename = bg_freq_filename[params[:bg]]
    session[:bgfreq_filename] = bgfreq_filename
    session[:bg] = params[:bg]

    if !bgfreq_filename.nil? && params[:bg] != "input_sequence"
      job.add_file(path(bg_freq_local_filename[params[:bg]]), working_directory)
      bgfreq_opt = " -B #{bgfreq_filename} "
    end

    job.command = <<-COMMAND.gsub(/\n/,'').squeeze(" ")
      java -jar #{working_directory}/plaac.jar -i #{filename} -c #{core_len} -a #{alpha} #{bgfreq_opt} > #{output_filename};
      touch #{working_directory}/results_ready
    COMMAND

    begin
      job.submit()
    rescue Exception => e
      @@log.error(e)
    end

    # If we only have one gene, skip to visualization
    if gene_count == 1
      params['token'] = job.token
      cookies[:prion_finder_picklist] = "1"
      display_prion_candidates()
    else
      redirect "#{$config[:apppath]}/candidates/#{job.token}"
    end
  end

  get '/candidates/:token/csv' do
    @job = Job.first(:token => params[:token])
    send_file File.join(@job.working_directory,"prion_candidates.tsv"),
      :type => 'text/plain',
      :filename => "prion_candidates.csv"
  end

  get '/candidates/:token' do
    @@log.info "candidates #{params[:token]}"
    @job = Job.first(:token => params[:token])
    # TODO: check that @job existed
    results_ready_file = File.join(@job.working_directory, "results_ready")
    prion_candidates_file = File.join(@job.working_directory, @@prion_candidates)


    # get header, and produce a hash of column symbol to array index.
    if File.exist?(results_ready_file)
      @@log.info "results ready"
      lines = File.open(prion_candidates_file).lines
      lines = lines.reject{|line| line =~ /^##/}
      header_line = lines.first
      if header_line.nil?
        redirect "#{$config[:apppath]}/candidates/#{@job.token}?bounce"
        return
      end
      @header = header_line.split(/\t/)
      @idx = Hash[@header.each_with_index.map{|column, i| [column.to_sym, i]}]

      # TODO: REVIEW: We should move this sorting into the java.
      sorted_flag_file = File.join(@job.working_directory, ".candiates.sorted")

      unless File.exist?(sorted_flag_file)
        i = 10
        until File.writable?(prion_candidates_file) && File.exist?(prion_candidates_file) || i == 0
          i -= 1
          sleep 1
        end
        @@log.info "sorting #{prion_candidates_file}"
        @output = lines.map{|line| line.split(/\t/)}
        # sort the candidates in reverse order from highest COREscore to lowest
        @output = @output.sort_by{|line| -1*line[@idx[:COREscore]].to_f rescue 0 }

        # save re-sorted file to make sure they are in the right order for next step
        File.open(prion_candidates_file, 'w') do |s|
          s.puts @header.join("\t")
          @output.each do |candidate|
            s.puts candidate.join("\t")
          end
        end

        FileUtils.touch sorted_flag_file
      end

      @candidates, @picklist = load_candidates(prion_candidates_file)
    else
      redirect "#{$config[:apppath]}/candidates/#{@job.token}"
      return
    end

    @@log.info "render!"
    haml :candidates
  end

  post '/visualize' do
    display_prion_candidates()
  end

  get '/visualize/:token/pdf' do
    @job = Job.first(:token => params[:token])
    send_file File.join(@job.working_directory,"prion_details.pdf"),
      :type => 'application/pdf',
      :filename => "prion_details.pdf"
  end

  get '/visualize/:token/strippdf' do
    @job = Job.first(:token => params[:token])
    send_file File.join(@job.working_directory,"prion_strip.pdf"),
      :type => 'application/pdf',
      :filename => "prion_strip.pdf"
  end

  get '/visualize/:token/strippng' do
    @job = Job.first(:token => params[:token])
    send_file File.join(@job.working_directory,"prion_strip.png"),
      :type => 'application/png',
      :filename => "prion_strip.png"
  end

  get '/visualize/:token/images/:id.png' do
    @job = Job.first(:token => params[:token])
    if @job.nil?
      @@log.warn("unknown image token #{params[:token]}")
      ""
    else
      image_filename="prion_details_#{"%05d" % [params[:id].to_i]}.png"
      file = File.join(@job.working_directory,image_filename)
      send_file file,
        :type => 'image/png',
        :filename => image_filename
    end
  end

  get '/visualize/:token/csv' do
    @job = Job.first(:token => params[:token])
    send_file File.join(@job.working_directory,@@prion_candidates_details),
      :type => 'text/plain',
      :filename => "prion_details.csv"
  end

  get '/visualize/:token' do
    @job = Job.first(:token => params[:token])
    candidates_filename = File.join(@job.working_directory,"prion_candidates_selected.tsv")
    @output = File.open(candidates_filename).readlines.map{|line| line.split(/\t/)}


    # TODO: this is stupid, we should be using the prion_finder_picklist to random access the file, not load the whole thing.
    prion_candidates_file = File.join(@job.working_directory, @@prion_candidates)
    header_line = File.open(prion_candidates_file).lines.first
    @header = header_line.try(:split, /\t/)
    @candidates = []
    unless @header.nil?
      @idx = Hash[@header.each_with_index.map{|column, i| [column.to_sym, i]}]
      all_candidates, @picklist = load_candidates(prion_candidates_file)
      selected_rows = cookies[:prion_finder_picklist].split(/,/).map{|i|i.to_i}
      all_candidates.each_with_index do |candidate, i|
        @candidates << candidate[0].split(/\t/) if selected_rows.include?(i)
      end
    end

    haml :visualize
  end

  #--------------------------------------------------------------------
  # Helpers

  # an absolute path for filename relative to @@base
  def path(filename)
    File.expand_path("#{filename}",@@base)
  end

  # writes filename with the contents from string
  def write_string(filename,string)
    File.open(filename,'wb+') do |f|
      f.write string
    end
  end

  # writes filename with the contents from iostream
  def write_file(filename,iostream)
    File.open(filename,'wb+') do |f|
      while (blk = iostream.read(65536))
          f.write blk
      end
    end
  end

  def time(comment='',&block)
    a = Time.now.to_f
    block.call
    b = Time.now.to_f
    p [comment,b-a]
  end

  def load_candidates(filename)
    @dataset = []
    @picklist = []

    page = params[:page] || 1
    filter = nil
    filter = /#{params[:filter]}/ unless [nil,''].include? params[:filter].nil?
    filtered = []
    picklist_ids = []
    cookie_picklist = cookies[:prion_finder_picklist]
    unless cookie_picklist.nil?
      begin
        cookie_picklist.gsub!(/[^0-9\,]/,'') # clean input
        picklist_ids = cookie_picklist.split(/,/).map(&:to_i)
      rescue Error => e
        puts "Could not load cookie list: #{e}"
        picklist_ids = []
      end
    end

    if filter
      # filter with grep
      time('ruby') {
        all_lines = File.open(filename).lines.each_with_index.map{|line,i| [line,i]}
        @picklist = all_lines.select{|pair| picklist_ids.include? pair[1] }
        filtered = all_lines.select{|pair| pair[0] =~ filter}
      }

      @total = filtered.size
    else
      @total = %x[ wc -l #{filename} ].strip.split(/ /).first.to_i
    end

    @per_page = 50

    @pages = (@total.to_f / @per_page).ceil

    @page = page.to_i || 1
    @page = 1 if @page > @pages

    start = (@page - 1) * @per_page
    start = 0 if start < 0

    if filter
      @dataset = filtered[start, @per_page]
    else
      # very fast if you want everything
      file = File.open(filename)
      all_lines = file.lines.each_with_index.map{|line,i| [line,i]}

      @picklist = all_lines.select{|pair| picklist_ids.include? pair[1] }

      (start).times { all_lines.next } # skip skip skip
      @dataset = all_lines.take(@per_page)
    end

    [@dataset, @picklist]
  end

  def add_gene_id(sequence)
    raw_sequences = []
    raw_sequence = ""
    sequence.strip.lines.each do |line|
      if line.strip == ""
        unless raw_sequence.empty?
          raw_sequences << raw_sequence
          raw_sequence = ""
        end
      else
        raw_sequence += line
      end
    end
    unless raw_sequence.empty?
      raw_sequences << raw_sequence
    end

    raw_sequences = raw_sequences.map{|sequence| sequence.strip}

    sequences = []
    raw_sequences.each do |sequence|
      if sequence[0] == '>'
        sequences << sequence
      else
        sequences << ">gene#{UUID.generate[0..7]}\n"+sequence
      end
    end
    sequences.join("\n\n")
  end

  def page_url(page)
    "?"+ { :filter => params[:filter], :page => page }.map{|e| e.join("=")}.join("&")
  end

  private

  def local(cmd)
    puts cmd
    p system(cmd)
  end

  def valid_fasta_file(filename)
    invalid_fasta_rex = /([^A-Z*-])/
    File.open(filename).each_with_index do |line,i|
      next if line =~ /^>/
      return [false, "Found invalid amino acid on line #{i}. '#{$1}'"] if line.strip =~ invalid_fasta_rex
    end
    return [true, nil]
  end

  def display_prion_candidates
    @@log.info "display_prion_candidates #{params.inspect}"
    @job = Job.first(:token => params["token"])
    @@log.debug @job
    @candidates = []

    prion_candidates_file = File.join(@job.working_directory, @@prion_candidates)

    core_len = session[:core_len]
    alpha = session[:alpha]

    # passed in from previous session
    bgfreq_filename = session[:bgfreq_filename]

    # rewrite prion_candidates, with selections
    session[:selected] = selected = [params[:selected]].flatten.compact.map(&:to_i)

    candidates_filename = File.join(@job.working_directory,"prion_candidates_selected.tsv")
    File.open(candidates_filename,'w') do |candidates_file|
      candidates, picklist = load_candidates(prion_candidates_file)
      picklist.each do |columns,i|
        name = columns.split(/\t/)[0]
        # TODO: ultimately the second "name" should be replaced by the a user-specified "common" gene name
        candidates_file.puts name + "\t" + name
        @candidates << name
      end
    end
    local("chmod ug+rw #{candidates_filename}")

    input_fasta = @@input_fasta
    output_filename = @@prion_candidates_details

    job = Job.new(@job.token)
    working_directory = job.working_directory

    # generate details
    if !bgfreq_filename.nil? && session[:bg] != "input_sequence"
      bg_freq_filename = ->(id){ id.nil? ? nil : "bg_freqs_#{id}.txt" }
      bg_freq_local_filename = ->(id){ id.nil? ? nil : "bg_freqs/#{bg_freq_filename[id]}" }
      job.add_file(path(bg_freq_local_filename[session[:bg]]), working_directory)
      bgfreq_opt = " -B #{bgfreq_filename} "
    end

    job.command = <<-COMMAND
    sh -c "cd #{@job.working_directory} &&
      java -jar ./plaac.jar -i #{input_fasta} -c #{core_len} -a #{alpha} -p #{candidates_filename} #{bgfreq_opt} > #{output_filename} &&
      ./plaac_plot.r prion_candidates_details.tsv prion_details.pdf &&
      ./plaac_plot.r prion_candidates_details.tsv prion_details.png &&
      ./plaac_plot.r prion_candidates_details.tsv prion_strip.pdf -c &&
      ./plaac_plot.r prion_candidates_details.tsv prion_strip.png -c ;
      touch details_ready
    "
    COMMAND
    #job.working_directory = @job.working_directory
    job.owner = session[:login]
    if job.save
      @@log.info "ok"
    else
      job.errors.each do |error|
        @@log.error(error)
      end
    end

    job.submit()

    redirect "#{$config[:apppath]}/visualize/#{job.token}"
  end

def read_bg_freqs
  bg_freqs = {}
  skipfirst = true
  File.open(@@background_meta).each do |line|
    if skipfirst
      skipfirst = false
      next
    end
    (uniprot_version, id_name, species_name, fasta_url, download_date, aa_freq_filename) = line.split(/\t+/)
    bg_freqs[id_name] = species_name
  end
  bg_freqs
end

  run! if __FILE__ == $0
end
